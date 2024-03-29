module functions
using Distributions
using LinearAlgebra
using DataStructures
using DelimitedFiles
using Random
using JuMP
using Gurobi
using Combinatorics
using DataFrames
using Dates
using DataFramesMeta
using BlackBoxOptim
# For monte carlo simulations, merger_cost_type = "num_of_targets_divided_by_buyer_size"

export expand_grid, shipmerger_struct, carrier_mod,
       #ship_merger_monte_carlo
       TB_toy,IS_outdomain,TB,IS_Sell,IS_Buyer,gen_X_mat,
	   gen_X_interaction, gen_merger_cost,
       gen_subsidy, gen_utility_matrix, solve_equilibrium, gen_data,
       score_b, maxscore_mc,
	   ship_merger_score_estimation,
	   gen_utility_est_without_subsidy,
	   gen_unmatched_utility_est,
	   gen_utility_est,
	   gen_subsidy_indicator_for_buyer_ii,
	   gen_merger_cost_for_buyer_ii,
	   extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk,
	   extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk,
	   extract_target_covariates_if_ii_is_buyer,
	   extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk,
	   extract_buyer_covariates_if_ii_is_buyer,
	   # score_b_est_data,
	   # score_bthis_full_X,
	   # score_bthis_scale_X_only,
	   # score_bthis_scale_and_scope_X_only,
	   # score_bthis_x45678_merger_cost
	   inner_construct_critical_value_CR_two_variables,
	   outer_construct_critical_value_CR_two_variables



function expand_grid(args...)
    nargs= length(args)
    if nargs == 0
      error("expand_grid need at least one argument")
    end
    iArgs= 1:nargs
    nmc= "Var" .* string.(iArgs)
    nm= nmc
    d= map(length, args)
    orep= prod(d)
    rep_fac= [1]
    # cargs = []
    if orep == 0
        error("One or more argument(s) have a length of 0")
    end
    cargs= Array{Any}(undef,orep,nargs)
    for i in iArgs
        x= args[i]
        nx= length(x)
        orep= Int(orep/nx)
        mapped_nx= vcat(map((x,y) -> repeat([x],y), collect(1:nx), repeat(rep_fac,nx))...)
        cargs[:,i] .= x[repeat(mapped_nx,orep)]
        rep_fac= rep_fac * nx
    end
    #convert(DataFrame,cargs)
	DataFrame(Tables.table(cargs, header = Symbol.(:x, axes(cargs, 2))))
end

mutable struct shipmerger_struct
    N::Int64        # number of carrier types in each market
    PN::Int64       # number of possible acquisition bundles in each marke
    X_mat::Array{Float64,2}  # realization of random vector [X_1, X_2, ..., X_S ]
	ton::Array{Float64,2} # tonnage size of each firm
    ϵ_mat::Array{Float64,2} # realization of random demand shocks: matrix size is N × ( N+1 )
	Bundle::Vector{String}  # auxiliary index objects    indx is 1,2,...,32 ⇒ bundle is vector of 0/1 say [0,0,0,1,1]
    eta::Float64 # measure of carriers of each type, here fixed to 1.
    β_0::Float64 # coefficient of covariates
	δ_0::Float64 # coefficient of subsidy
	γ_0::Float64 # coefficient of additional mergers
	ton_dim::Int64 # dimension of carrier type for each firm
end

function carrier_mod(;N::Int64=4,#11,#10,#,9,
                    ton_dist::UnivariateDistribution=LogNormal(2.0,1.0),  # distribution for size
                    ϵ_dist::UnivariateDistribution=Normal(0.0,1.0),  # distribution for demand shocks
                    eta::Float64=1.0,
                    β_0::Float64=1.0, #0.0
					δ_0::Float64=2.0,
					γ_0::Float64=2.0,
                    randomseed::Int64 = 1,
					ton_dim::Int64 = 2
                    )
    Random.seed!(randomseed)
    PN = 2^N;
    eps_dims = N + 1;   # dimension of ϵ for each individual
    #ton_mat = rand(ton_dist,N,ton_dim)
	ton_mat = rand(ton_dist,N,ton_dim)./100
	X_mat = ton_mat
    #ω_vec = rand(ω_dist,N)
    ϵ_post_mat = rand(ϵ_dist,N,PN) # efficient version
    # generate auxiliary Bundle vector to indx acquire target
    Bundle = Array{String}(undef,PN)
    for i =1:PN
        Bundle[i] = string(i-1,base=2,pad=N)
    end
    m = shipmerger_struct(N,PN,X_mat,ton_mat,ϵ_post_mat,Bundle,eta,β_0,δ_0,γ_0,ton_dim)
    return m
end

function TB_toy(N::Integer,w::Integer)  # N is number of players, but w is id of bundle.
    PN = 2^N;
    Bundle = Array{String}(undef, PN)
    for i = 1:PN
        Bundle[i] = string(i-1,base=2,pad=N)
    end
    cache_TB = Vector{Int64}(undef,N)
    for i=1:N
        cache_TB[i] = parse(Int,(Bundle[w])[i])
    end
    return cache_TB
end

# indicator to hard code utility to reach integer equilibrium
function IS_outdomain(N::Integer,i::Integer,ac::Integer)
    if sum(TB_toy(N,ac)) != 1
        if TB_toy(N,ac)[i] == 1
            return true
        else
            return false
        end
    else
        return false
    end
end

function TB(m::shipmerger_struct, w::Integer)    # (target indx) -> target bundle
    N = m.N
    BD = m.Bundle
    cache_TB = Vector{Int64}(undef,N)
    for i=1:N
        cache_TB[i] = parse(Int,(BD[w])[i])
    end
    return cache_TB
end
# function that return true if t choose ac means t sell itself; return false if t acquire ac.
function IS_Sell(m::shipmerger_struct,t::Integer,ac::Integer)
    if sum(TB(m,ac)) == 1
        if TB(m,ac)[t] ==1
            return true
        else
            return false
        end
    else
        return false
    end
end

function IS_Buyer(m::shipmerger_struct, t::Integer, ac::Integer)
	if IS_Sell(m,t,ac) == 1
        return false # t is a seller
    else
		if TB(m,ac) == zeros(m.N)
			return false # t is unmatched
		else
			return true # t is a buyer
		end
    end
end

function gen_X_mat(m::shipmerger_struct)
	ton = m.ton
	ton_dim = m.ton_dim
	X_dim = ton_dim*2 # size and share
	sum_ton = zeros(m.N, 1)
	for i = 1:m.N
		sum_ton[i] = sum(ton[i,:])
	end
	X_mat = hcat(sum_ton,ton)
	#X_mat = hcat(sum_ton,ton)/.100
	X_share_mat = X_mat[:,2:end]./X_mat[:,1]
	# 1st column is total tonnage, 2nd column gives log-size, and 3rd column gives share
	X_mat = hcat(X_mat, X_share_mat)
	return X_mat
end

function gen_X_interaction(m::shipmerger_struct, X_mat::Array{Float64,2})
	dim_X_mat = size(X_mat)[2]
	dim_X_mat_share = Int((dim_X_mat - 1)/2)
	X_interaction = zeros(m.N, m.PN, dim_X_mat)
	target_X = zeros(m.PN, dim_X_mat)
	buyer_X = copy(X_mat)
	for j = 1:m.PN
		temp_index = TB(m,j)
		for k = 1:m.N
			# total tonnage(column 1)
			target_X[j,1] = target_X[j,1] + temp_index[k].*buyer_X[k,1]
			# types of tonnage and logalization(column 2,3)
			for l = 2:1+dim_X_mat_share
			    target_X[j,l] = target_X[j,l] + temp_index[k].*buyer_X[k,l]
			end
		end
		target_X[j,2:1+dim_X_mat_share] # = log.(target_X[j,2:2+dim_X_mat_share])
		# types of shares
		enum = sum(temp_index.*buyer_X[:,2+dim_X_mat_share:end],dims=1)
		denom = sum(temp_index)
		target_X[j,2+dim_X_mat_share:end] = enum/denom
	end
	# taking log
	if 0 == 1
	    target_X[2:end,1:1+dim_X_mat_share] = copy.(log.(1 .+target_X[2:end,1:1+dim_X_mat_share]))
	    buyer_X[1:end,1:1+dim_X_mat_share] = copy.(log.(1 .+buyer_X[1:end,1:1+dim_X_mat_share]))
	end
	# add exogenous shocks to population matrix to avoid non-integer solution
	# pop_shock = rand(m.N,m.PN,m.ton_dim)
	for i = 1:m.N, j = 1:m.PN
		X_interaction[i,j,:] = buyer_X[i,:].*target_X[j,:]# + pop_shock[i,j,:]
	end
	return X_interaction, target_X
end
# For monte carlo simulations, merger_cost_type = "num_of_targets_divided_by_buyer_size"
function gen_merger_cost(m::shipmerger_struct;
	                     merger_cost_type = "only_num_of_targets")
	ton = m.ton
	buyer_X = gen_X_mat(m)
	merger_cost = zeros(m.N,m.PN)
	if merger_cost_type == "only_num_of_targets"
		for j = 1:m.PN
			temp_index = TB(m,j)
			for k = 1:m.N
				# merger cost
				num_of_firms_in_coalition = sum(temp_index)
				merger_cost[k,j] = num_of_firms_in_coalition
			end
		end
	elseif merger_cost_type == "num_of_targets_divided_by_buyer_size"
		for j = 1:m.PN
			temp_index = TB(m,j)
			for k = 1:m.N
				# merger cost
				num_of_firms_in_coalition = sum(temp_index)
				total_tonnage = sum(ton[k,:])
				merger_cost[k,j] = num_of_firms_in_coalition/log(total_tonnage*100)
			end
		end
	end
	return merger_cost
end

function gen_subsidy(m::shipmerger_struct, threshold_tonnage::Any, subsidy_amount::Any; subsidy_type::Any)
	ton = m.ton
	subsidy_index_mat = zeros(m.N,m.PN)
	total_tonnage_in_coalition = zeros(m.N,m.PN)
	for j = 1:m.PN
		temp_index = TB(m,j)
		total_tonnage_in_coalition[:,j] = sum(ton,dims = 2).+sum(temp_index.*ton)
		for i = 1:m.N
			if IS_Sell(m,i,j)
				total_tonnage_in_coalition[i,j] = sum(ton,dims = 2)[i]
			end
			# subsidy
			if total_tonnage_in_coalition[i,j] >= threshold_tonnage
				if subsidy_type == "to_buyer"
					subsidy_index_mat[i,j] = 1*subsidy_amount
				elseif subsidy_type == "shared"
					subsidy_index_mat[i,j] = subsidy_amount/sum(temp_index)
				end
			end
		end
	end
	return subsidy_index_mat, total_tonnage_in_coalition
end

function gen_utility_matrix(m::shipmerger_struct, X_interaction::Any, threshold_tonnage::Any, subsidy_amount::Any; subsidy_type::Any)
	β = m.β_0
	γ = m.γ_0
	δ = m.δ_0
	utility = zeros(m.N,m.PN)
	merger_cost = gen_merger_cost(m)
	subsidy_index_mat,total_tonnage_in_coalition = gen_subsidy(m,
	                                                           threshold_tonnage,
	                                                           subsidy_amount,
															   subsidy_type = subsidy_type)
	#exogenous_shock = randn(m.N, m.PN)
	for i = 1:m.N, j= 1:m.PN
		if IS_Sell(m,i,j)
			# rescale selling payoff from quadratic interaction (x_i*x_i) into (x_i*1)
			utility[i,j] = sqrt(X_interaction[i,j,2])*1 + sqrt(X_interaction[i,j,4])*β +
			               subsidy_index_mat[i,j]*δ + m.ϵ_mat[i,j]
		else
			utility[i,j] = X_interaction[i,j,2]*1 + X_interaction[i,j,4]*β +
			               subsidy_index_mat[i,j]*δ + m.ϵ_mat[i,j]
		end
        utility[i,j] = utility[i,j] #+ exogenous_shock[i,j]
	end
	# assigning negative infinity to unreal allocation. (ex.) (1, (101)) and (1,(1,1,1))
	for i = 1:m.N
        for j = 1:m.PN
            if IS_outdomain(m.N,i,j) # hard code an continuum agent model to reach an integer equilibrium
                utility[i,j] = -99999999 # unreal payoff
            elseif IS_Sell(m,i,j)
                utility[i,1] = copy(utility[i,j]) # unmatched payoff(=buy myself now)
				utility[i,j] = 0 + m.ϵ_mat[i,j]  # selling payoff
            else
				utility[i,j] = utility[i,j] - merger_cost[i,j]*γ
			end
        end
    end
	return utility
end

function solve_equilibrium(m::shipmerger_struct,threshold_tonnage::Any,subsidy_amount::Any; subsidy_type::Any)
	X_interaction,target_X = gen_X_interaction(m, gen_X_mat(m))
	utility = gen_utility_matrix(m,
	                             X_interaction,
	                             threshold_tonnage,
								 subsidy_amount,
								 subsidy_type = subsidy_type)
	N = m.N
	PN = m.PN
	eta = m.eta
	#env = Gurobi.GRBenv()
    model = JuMP.Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", 100)
    set_optimizer_attribute(model, "Presolve", 0)
    JuMP.@variable(model, 0<=x[i=1:N,j=1:PN]<=1)
	#JuMP.@variable(model, x[i=1:N,j=1:PN], binary=true)
	# define a function to write the expression for buyers of type i=1,...,N. this expr will appear in double coinsidence condition.
	function buyer_expression_ipopt(m::shipmerger_struct, t::Integer)
		isthis(x,y) = x == y
		sum_x_expression = 0
		for i=1:m.N
			for j=1:m.PN
				if IS_Sell(m,i,j)
					nothing
				else
					ac_target = findall(isone,TB(m,j))
					ac_bool = sum([isthis(t,k) for k in ac_target])
					if ac_bool>=1 # it means x[m,n] is buyer of type t
						sum_x_expression = sum_x_expression + x[i,j]
					else
						nothing
					end
				end
			end
		end
		return sum_x_expression
	end
	println("# feasibility constraints(double-coinside const):\n")
	mm = m
	for m = 1:N
		for n = 1:PN
			if IS_Sell(mm,m,n)
				#double_cc_expression_0 ="subject to feas_mu_$m:"
				double_cc_expression_1 = x[m,n]# =
				double_cc_expression_2 = buyer_expression_ipopt(mm,m)
				#double_cc = double_cc_expression_1-double_cc_expression_2
				#write(f, double_cc)
				@constraint(model,double_cc_expression_1 == double_cc_expression_2)
			else
				nothing
			end
		end
		#@show m
	end
	#println( "# feasibility constraints(measure const):\n")
	@constraint(model, feas_rho[i=1:N], sum(x[i,j] for j in 1:PN)== eta)

    JuMP.@objective(model, Max, sum(x[i,j]*utility[i,j] for i in 1:N, j in 1:PN))
    #println("Time for optimizing model:")
    @time JuMP.optimize!(model)
    # show results
    objv = JuMP.objective_value(model)
    #println("objvalue　= ", objv)
    #matches = round(mylp$solution, 1)
    matches = JuMP.value.(x)
	matches = round.(matches, digits = 1)
    return utility, matches
end

function gen_data(m::shipmerger_struct;threshold_tonnage = 100,subsidy_amount = 100, subsidy_type = "to_buyer")
	utility, matches = solve_equilibrium(m, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
	round.(matches, digits = 2)
	buyer_X = gen_X_mat(m)
	X_interaction, target_X = gen_X_interaction(m, buyer_X)
	buydata = buyer_X
	buyid = Array{Int64,1}(1:m.N)
	buydata = hcat(buyid, buydata)
	#buydata = convert(DataFrame, buydata)
	buydata = DataFrame(Tables.table(buydata, header = Symbol.(:x, axes(buydata, 2))))
	x = ["size"]
	y = ["share"]
	z = ["total_size_buyer"]
	for i = 1:m.ton_dim
		z = vcat(z, )
	end
	for i = 1:m.ton_dim
		z = vcat(z, map(*, x, "$(i)"))
	end
	for i = 1:m.ton_dim
		z = vcat(z, map(*, y, "$(i)"))
	end
	rename!(buydata, [:id, :total_size_buyer, :size1_buyer, :size2_buyer, :share1_buyer, :share2_buyer])
	#tardata = mvrnorm(I, mu=means, Sigma=covars)
	tardata = target_X
	tarid = Array(1:m.PN)
	tardata = hcat(tarid, tardata)
	#tardata = convert(DataFrame, tardata)
	tardata = DataFrame(Tables.table(tardata, header = Symbol.(:x, axes(tardata, 2))))
	rename!(tardata, [:id, :total_size_target, :size1_target, :size2_target, :share1_target, :share2_target])
	matchmaker = expand_grid(buyid, tarid)
	rename!(matchmaker, [:buyid, :tarid])
	matchdat = DataFrames.leftjoin(matchmaker, tardata,
                               on = [:tarid => :id])
	matchdat = DataFrames.leftjoin(matchdat, buydata,
	                           on = [:buyid => :id])
	sort!(matchdat, [:buyid, :tarid]);
	#matchdat = within(matchdat, mval <- mval + rnorm(length(matchdat$mval), mean = 0, sd_err) )
	mval = vec(utility') .+ vec(m.ϵ_mat')
	matchdat = hcat(matchdat, mval)
	rename!(matchdat, :x1 => :mval)
	matchdat = hcat(matchdat, vec(matches'))
	rename!(matchdat, :x1 => :matches)
	DataFramesMeta.@linq obsd = matchdat |>
	  where(:matches .> 0.0)
	  return utility, obsd
end

function score_b(m::shipmerger_struct, beta::Vector{Float64},
	             data::DataFrame, num_agents::Int64,
				 threshold_tonnage::Any, subsidy_amount::Any;
				 subsidy_type::Any,
				 compare_with_and_without_subsidy = "no")
    I = num_agents
    temp = [Combinatorics.combinations(1:I,2)...]
    index_list = Array{Int64,2}(undef, length(temp), 2)
    for i in 1:length(temp)
        index_list[i,1] = temp[i][1]
        index_list[i,2] = temp[i][2]
    end
	dataX = convert(Array{Float64,2},data[:,2+5+1:2+5+5]) # pick up buyer covariates
	X_interaction_data, target_X = gen_X_interaction(m, dataX)
	merger_cost = gen_merger_cost(m)
	subsidy_index_mat, total_tonnage_in_coalition = gen_subsidy(m, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
	#X_interaction[i,j,2]*1 + X_interaction[i,j,4]*β + subsidy_index_mat[i,j]*subsidy_amount*δ - merger_cost[i,j]*γ + m.ϵ_mat[i,j]
	beta_hat, delta_hat, gamma_hat = beta # for Optim
	utility = zeros(m.N, m.PN)
	utility_without_subsidy = zeros(m.N, m.PN)
	#exogenous_shock = randn(m.N, m.PN)
	#println("******
	#Assign utility specification
	#*******")
	for i = 1:m.N, j= 1:m.PN
		utility[i,j] = X_interaction_data[i,j,2]*1.0 + X_interaction_data[i,j,4]*beta_hat + subsidy_index_mat[i,j]*delta_hat + m.ϵ_mat[i,j]
        #utility[i,j] = utility[i,j] + exogenous_shock[i,j]
	end
	# assigning negative infinity to unreal allocation. (ex.) (1, (101)) and (1,(1,1,1))
	for i = 1:m.N
        for j = 1:m.PN
            if IS_outdomain(m.N,i,j)      # hard code an continuum agent model to reach an integer equilibrium
                utility[i,j] = -9999999
			elseif IS_Sell(m,i,j)
                utility[i,1] = copy(utility[i,j]) # unmatched payoff
				utility[i,j] = m.ϵ_mat[i,j]  # selling payoff
            else
				utility[i,j] = utility[i,j] - merger_cost[i,j]*gamma_hat
			end
        end
    end
	for i = 1:m.N, j= 1:m.PN
		utility_without_subsidy[i,j] = copy(utility[i,j])
		utility_without_subsidy[i,j] = utility_without_subsidy[i,j] - subsidy_index_mat[i,j]*delta_hat
	end
	# construct inequality
	matching_index = hcat(data.buyid, data.tarid)
	unmatched_vec = zeros(m.N)
	ineq = zeros(length(index_list[:,1]), 20, 20)
    for i = 1:length(index_list[:,1])
		idx = index_list[i,:] # pick buyer id pair
		target_bundle_id = [matching_index[idx[1],2], matching_index[idx[2],2]]
		if IS_Buyer(m,idx[1],target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
			#println("iter $i = Case 1: both firms are buyers.")
			#First, I construct matching maximum score inequality for two acquiring firms without price data:
			payoff_obs_match1 = utility[idx[1],target_bundle_id[1]] # buyer
			payoff_obs_match2 = utility[idx[2],target_bundle_id[2]] # buyer
			k_index = zeros(Int64,0)
			for k = 1:m.N # extract member of firm 1's coalition
				m.Bundle[target_bundle_id[1]][k]== '1'
				if m.Bundle[target_bundle_id[1]][k]== '1'
				    k_index = vcat(k_index, k)
				end
			end
			h_index = zeros(Int64,0)
			for h = 1:m.N # extract member of firm 2's coalition
				m.Bundle[target_bundle_id[2]][h]== '1'
				if m.Bundle[target_bundle_id[2]][h]== '1'
				    h_index = vcat(h_index, h)
				end
			end
			for kk = 1:length(k_index), hh = 1:length(h_index)
				k = k_index[kk]
				h = h_index[hh]
				swapped_bundle_id1 = m.Bundle[target_bundle_id[1]][1:k-1]*repeat("0", 1)*m.Bundle[target_bundle_id[1]][k+1:m.N] # drop firm k
				swapped_bundle_id1 = swapped_bundle_id1[1:h-1]*repeat("1", 1)*swapped_bundle_id1[h+1:m.N] # add firm h
				swapped_bundle_id2 = m.Bundle[target_bundle_id[2]][1:h-1]*repeat("0", 1)*m.Bundle[target_bundle_id[2]][h+1:m.N] # drop firm h
				swapped_bundle_id2 = swapped_bundle_id2[1:k-1]*repeat("0", 1)*swapped_bundle_id2[k+1:m.N] # add firm k
				payoff_unobs_match1 = utility[idx[1],Int.(findall(m.Bundle .== swapped_bundle_id1))[1]] # drop firm k
				payoff_unobs_match2 = utility[k,Int.(findall(m.Bundle .== swapped_bundle_id2))[1]]
				ineq[i, kk, hh] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			end
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				if subsidy_index_mat[idx[1],target_bundle_id[1]] == 1
					# the observed pairwise data satisfy the subsidy threshold
				    ineq[i, 1, length(h_index) + 1] = utility[idx[1],1] - utility_without_subsidy[idx[1],target_bundle_id[1]] # without subsidy
					ineq[i, 1, length(h_index) + 2] = payoff_obs_match1 - utility[idx[1],1] # with subsidy
				end
				if subsidy_index_mat[idx[2],target_bundle_id[2]] == 1
					# the observed pairwise data satisfy the subsidy threshold
					ineq[i, 1, length(h_index) + 3] = utility[idx[2],1] - utility_without_subsidy[idx[2],target_bundle_id[2]]
					ineq[i, 1, length(h_index) + 4] = payoff_obs_match2 - utility[idx[2],1] # with subsidy
				end
			end
		elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && IS_Sell(m,idx[2],target_bundle_id[2])
			#println("iter $i = Case 2: firm 1 is a buyer and firm 2 is a seller.")
			#Second, I construct inequalities from an observed coalition:
			payoff_obs_match1 = utility[idx[1],target_bundle_id[1]] # buyer
			# choose a firm out of coalition of buyer a
			k_index = zeros(Int64,0)
			for k = 1:m.N
				m.Bundle[target_bundle_id[1]][k]== '1'
				if m.Bundle[target_bundle_id[1]][k]== '1'
				    k_index = vcat(k_index, k)
				end
			end
			for kk = 1:length(k_index)
				k = k_index[kk]
				swapped_bundle_id1 = m.Bundle[target_bundle_id[1]][1:k-1]*repeat("0", 1)*m.Bundle[target_bundle_id[1]][k+1:m.N] # drop firm k
				#swapped_bundle_id2 = repeat("0", m.N) # unmatched firm k
				payoff_unobs_match1 = utility[idx[1],Int.(findall(m.Bundle .== swapped_bundle_id1))[1]] # drop firm k
				#payoff_unobs_match2 = utility[k,Int.(findall(m.Bundle .== swapped_bundle_id2))[1]] # unmatched firm k
				ineq[i,kk,1] = payoff_obs_match1 - payoff_unobs_match1 #- payoff_unobs_match2
			end
			if compare_with_and_without_subsidy == "yes"
			    # compare utility with and without subsidy
			    if subsidy_index_mat[idx[1],target_bundle_id[1]] == 1
				# the observed pairwise data satisfy the subsidy threshold
			    ineq[i, 1, 2] = utility[idx[1],1] - utility_without_subsidy[idx[1],target_bundle_id[1]] #without subsidy
				ineq[i, 1, 3] = payoff_obs_match1 - utility[idx[1],1] # with subsidy
			    end
		    end
		elseif IS_Sell(m,idx[1],target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
			#println("iter $i = Case 3: firm 1 is a seller and firm 2 is a buyer.")
			#Second, I construct inequalities from an observed coalition:
			payoff_obs_match2 = utility[idx[2],target_bundle_id[2]] # buyer
			# choose a firm from a coalition
			k_index = zeros(Int64,0)
			for k = 1:m.N
				#m.Bundle[target_bundle_id[2]][k]== '1'
				if m.Bundle[target_bundle_id[2]][k]== '1'
				    k_index = vcat(k_index, k)
				end
			end
			for kk = 1:length(k_index)
				k = k_index[kk] # pick up a firm which will be apart from a coaition
				#swapped_bundle_id1 = repeat("0", m.N) # unmatched firm k
				swapped_bundle_id2 = m.Bundle[target_bundle_id[2]][1:k-1]*repeat("0", 1)*m.Bundle[target_bundle_id[2]][k+1:m.N] # drop firm k
				#payoff_unobs_match1 = utility[k,Int.(findall(m.Bundle .== swapped_bundle_id1))[1]] # unmatched firm k
				payoff_unobs_match2 = utility[idx[2],Int.(findall(m.Bundle .== swapped_bundle_id2))[1]] # drop firm k
				ineq[i,kk,1] = payoff_obs_match2 - payoff_unobs_match2
			end
			if compare_with_and_without_subsidy == "yes"
			    # compare utility with and without subsidy
			    if subsidy_index_mat[idx[2],target_bundle_id[2]] == 1
				    # the observed pairwise data satisfy the subsidy threshold
			        ineq[i, 1, 2] = utility[idx[2],1] - utility_without_subsidy[idx[2],target_bundle_id[2]]
					ineq[i, 1, 3] = payoff_obs_match2 - utility[idx[2],1]
			    end
			end
		elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
			#println("iter $i = Case 4: firm 1 is a buyer and firm 2 is unmatched.")
			#Third, I construct inequalities from an unmatched target:
			payoff_obs_match1 = utility[idx[1],target_bundle_id[1]] # buyer
			payoff_obs_match2 = utility[idx[2],target_bundle_id[2]] # unmatched
			# choose a firm out of coalition
			k_index = zeros(Int64,0)
			for k = 1:m.N
				m.Bundle[target_bundle_id[1]][k]== '1'
				if m.Bundle[target_bundle_id[1]][k]== '1'
				    k_index = vcat(k_index, k)
				end
			end
			for kk = 1:length(k_index)
				k = k_index[kk]
				b = idx[2]
				swapped_bundle_id1 = m.Bundle[target_bundle_id[1]][1:k-1]*repeat("0", 1)*m.Bundle[target_bundle_id[1]][k+1:m.N] # drop firm k
				swapped_bundle_id1 = swapped_bundle_id1[1:b-1]*repeat("1", 1)*swapped_bundle_id1[b+1:m.N] # add firm b
				swapped_bundle_id2 = repeat("0", b-1)*repeat("1", 1)*repeat("0", m.N-b) # sell firm b
				payoff_unobs_match1 = utility[idx[1],Int.(findall(m.Bundle .== swapped_bundle_id1))[1]] # drop firm k and add firm b
				payoff_unobs_match2 = utility[b,Int.(findall(m.Bundle .== swapped_bundle_id2))[1]] # sell firm b
				ineq[i,kk,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			end
			if compare_with_and_without_subsidy == "yes"
			    # compare utility with and without subsidy
			    if subsidy_index_mat[idx[1],target_bundle_id[1]] == 1
				    # the observed pairwise data satisfy the subsidy threshold
			        ineq[i, 1, 2] = utility[idx[1],1] - utility_without_subsidy[idx[1],target_bundle_id[1]]
					ineq[i, 1, 3] = payoff_obs_match1 - utility[idx[1],1]
			    end
			end
		elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
			#println("iter $i = Case 5: firm 1 is unmatched and firm 2 is a buyer.")
			#Third, I construct inequalities from an unmatched target:
			payoff_obs_match1 = utility[idx[1],target_bundle_id[1]] # unmatched
			payoff_obs_match2 = utility[idx[2],target_bundle_id[2]] # buyer
			# choose a firm out of coalition
			k_index = zeros(Int64,0)
			for k = 1:m.N
				m.Bundle[target_bundle_id[2]][k]== '1'
				if m.Bundle[target_bundle_id[2]][k]== '1'
				    k_index = vcat(k_index, k)
				end
			end
			for kk = 1:length(k_index)
				k = k_index[kk]
				b = idx[1]
				swapped_bundle_id1 = repeat("0", b-1)*repeat("1", 1)*repeat("0", m.N-b) # sell firm b
				swapped_bundle_id2 = m.Bundle[target_bundle_id[2]][1:k-1]*repeat("0", 1)*m.Bundle[target_bundle_id[2]][k+1:m.N] # drop firm k
				swapped_bundle_id2 = swapped_bundle_id2[1:b-1]*repeat("1", 1)*swapped_bundle_id2[b+1:m.N] # add firm b
				payoff_unobs_match1 = utility[b,Int.(findall(m.Bundle .== swapped_bundle_id1))[1]] # sell firm b
				payoff_unobs_match2 = utility[idx[2],Int.(findall(m.Bundle .== swapped_bundle_id2))[1]] # drop firm k and add firm b
				ineq[i,kk,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			end
			#if compare_with_and_without_subsidy == "yes"
			    # compare utility with and without subsidy
			#    if subsidy_index_mat[idx[2],target_bundle_id[2]] == 1
			#	    # the observed pairwise data satisfy the subsidy threshold
			#        ineq[i, 1, 2] = utility[idx[2],1] - utility_without_subsidy[idx[2],target_bundle_id[2]]
			#    end
			#end
	    elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
			#println("iter $i = Case 6: both picked firms are unmatched.")
			# Fourth, I construct inequalities from IR conditions without subsidy
			# Here, deviation comes from merging unmatched pairwise firm
			payoff_obs_match1 = utility[idx[1],target_bundle_id[1]] # unmatched
			payoff_obs_match2 = utility[idx[2],target_bundle_id[2]] # unmatched
			swapped_bundle_id1 = repeat("0", idx[2]-1)*repeat("1", 1)*repeat("0", m.N-idx[2]) # match with unmatched firm 2
			swapped_bundle_id2 = repeat("0", idx[1]-1)*repeat("1", 1)*repeat("0", m.N-idx[1]) # match with unmatched firm 1
			payoff_unobs_match1 = utility[idx[1],Int.(findall(m.Bundle .== swapped_bundle_id1))[1]] # swapped
			payoff_unobs_match2 = utility[idx[2],Int.(findall(m.Bundle .== swapped_bundle_id2))[1]] # swapped
			ineq[i,1,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
		else
			#println("iter $i = Case 7,8,9: no firms are buyers.")
			ineq[i,1,1] = 0
		end
    end
	total_num_ineq = sum(ineq .!= 0)
    res = sum(ineq.>0)
    return res, total_num_ineq, utility
end

function maxscore_mc(num_agents::Int64;
                     num_its::Int64 = 10,
                     threshold_tonnage = 10,
					 subsidy_amount = 10,
					 subsidy_type = "to_buyer")
    start = Dates.unix2datetime(time())
	param_dim = 3
    myests = zeros(num_its, param_dim)
	mynum_correct_ineq = zeros(num_its, 1)
	truenum_correct_ineq = zeros(num_its, 1)
	num_all_ineq = zeros(num_its)
	solved_case_index = zeros(Int64, 0)
    for i = 1:num_its
       println("Create obsdat for iteration $i \n" )
	   @time m = carrier_mod(randomseed = i,
							N = num_agents,
							β_0 = 1.0,# coefficient of covariates
	  					    δ_0 = 1.0,# coefficient of subsidy
	  					    γ_0 = 2.0,# coefficient of additional merger cost
							ton_dim= 2)
	   #matches = solve_equilibrium(m,gen_utility_matrix)
	   #utility, matches = solve_equilibrium(m, threshold_tonnage)
	   utility, obsd = gen_data(m, threshold_tonnage = threshold_tonnage, subsidy_amount = subsidy_amount, subsidy_type = subsidy_type)
	   matching_index = hcat(obsd.buyid, obsd.tarid)
	   if size(obsd)[1] > m.N
		   	#println("skip because the matching outcome is not integer")
	   else
		   solved_case_index = vcat(solved_case_index, i)
		   function score_bthis(beta::Vector{Float64}, obsdat::DataFrame)
			   score, ~, ~ = score_b(m, beta, obsdat, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
			   res = -1.0*score + 100000.0 # need to be Float64 for bboptimize
			   #println("score:$(res) \n" )
			   return res
		   end
		   m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsd); SearchRange = (0.0, 10.0), NumDimensions = param_dim, Method = :de_rand_1_bin,  MaxSteps=10)
		   #println("score: ", score_bthis(m_res.archive_output.best_candidate, obsd))
		   mynum_correct_ineq[i,1],~ = score_b(m, m_res.archive_output.best_candidate, obsd, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
		   truenum_correct_ineq[i,1],total_num_ineq = score_b(m, [m.β_0, m.δ_0, m.γ_0], obsd, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
		   #println("the number of correct inequalities: ", mynum_correct_ineq[i,1])
		   #println("the TRUE number of score of correct: ", truenum_correct_ineq[i,1])
		   myests[i,:] = m_res.archive_output.best_candidate
		   num_all_ineq[i] = total_num_ineq
	   end
   end
   #runtime = proc.time() - ptm
   finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
   #mean = mean(vec-meanmat)
   #sqrt = sqrt(mean((vec-meanmat)^2))
   meanmat = reshape(repeat([m.β_0, m.δ_0, m.γ_0], num_its), param_dim, num_its)' # true parameter
   res_abs_mean_err = mean(abs.((myests.-meanmat)[solved_case_index,:]), dims = 1)
   res_sqrt = mean(((myests.-meanmat)[solved_case_index,:]).^2, dims = 1)
   mynum_correct_ineq = mynum_correct_ineq[solved_case_index,:]
   truenum_correct_ineq = truenum_correct_ineq[solved_case_index,:]
   num_all_ineq = num_all_ineq[solved_case_index]
   #@show myests
   return myests, res_abs_mean_err, res_sqrt, mynum_correct_ineq, truenum_correct_ineq, num_all_ineq
end

#------------------------------------#
# used in ship_merger_score_estimation
#------------------------------------#

function extract_buyer_covariates_if_ii_is_buyer(ii, data)
	# pick individual own covariates
	buyer_X_scale = Vector(data[ii,5:8])
	buyer_X_scope = Vector(data[ii,11:14])
	res = vcat(buyer_X_scale,buyer_X_scope)
	return res
end

function extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(ii, kk, data)
	# pick individual own covariates
	# pick up group names
	group_names = data[ii,4]
	# pick up deviating member
	subdata = @linq data |>
	    where(:group .== group_names)
	deviater_X_scale = Vector(subdata[kk,5:8])
	deviater_X_scope = Vector(subdata[kk,11:14])
	res = vcat(deviater_X_scale, deviater_X_scope)
	return res
end

function extract_target_covariates_if_ii_is_buyer(ii, data; info_sum = temp_info_sum)
	#global liner_sum, special_sum, tanker_sum, tramper_sum
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
	# pick up group names
	group_names = data[ii,4]
	# pick individual own covariates
	buyer_X_scale = Vector(data[ii,5:8])
	buyer_total_sum = data.total[ii]
	target_total_sum = data.group_total_tonnage[ii]
	# construct scale covariates
	target_liner_sum = zeros(0)
	target_special_sum = zeros(0)
	target_tramper_sum = zeros(0)
	target_tanker_sum = zeros(0)
	for iter = 1:size(liner_sum,1)
		if liner_sum.group[iter] == group_names
			target_liner_sum = liner_sum.liner_sum[iter]
		end
		if special_sum.group[iter] == group_names
			target_special_sum = special_sum.special_sum[iter]
		end
		if tramper_sum.group[iter] == group_names
			target_tramper_sum = tramper_sum.tramper_sum[iter]
		end
		if tanker_sum.group[iter] == group_names
			target_tanker_sum = tanker_sum.tanker_sum[iter]
		end
	end
	# construct scope covariates
	target_liner_share = (target_liner_sum-buyer_X_scale[1])/(target_total_sum-buyer_total_sum)
	target_special_share = (target_special_sum-buyer_X_scale[2])/(target_total_sum-buyer_total_sum)
	target_tramper_share = (target_tramper_sum-buyer_X_scale[3])/(target_total_sum-buyer_total_sum)
	target_tanker_share = (target_tanker_sum-buyer_X_scale[4])/(target_total_sum-buyer_total_sum)
	res = vcat(target_liner_sum, target_special_sum, target_tramper_sum, target_tanker_sum,
	       target_liner_share, target_special_share, target_tramper_share, target_tanker_share)
	return res
end

function extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(ii, kk, data; info_sum = temp_info_sum)
	#global liner_sum, special_sum, tanker_sum, tramper_sum
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
	# pick up group names
	group_names = data[ii,4]
	# pick up deviating member
	subdata = @linq data |>
	    where(:group .== group_names)
	subdata.id[kk] == data.id[ii]
	deviater_X_scale = Vector(subdata[kk,5:8])
	deviater_total_sum = subdata[kk,9]
	# pick individual own covariates
	buyer_X_scale = Vector(data[ii,5:8])
	buyer_total_sum = data.total[ii]
	target_total_sum = data.group_total_tonnage[ii]
	# construct scale covariates
	target_liner_sum = zeros(0)
	target_special_sum = zeros(0)
	target_tramper_sum = zeros(0)
	target_tanker_sum = zeros(0)
	for iter = 1:size(liner_sum,1)
		if liner_sum.group[iter] == group_names
			target_liner_sum = liner_sum.liner_sum[iter]
		end
		if special_sum.group[iter] == group_names
			target_special_sum = special_sum.special_sum[iter]
		end
		if tramper_sum.group[iter] == group_names
			target_tramper_sum = tramper_sum.tramper_sum[iter]
		end
		if tanker_sum.group[iter] == group_names
			target_tanker_sum = tanker_sum.tanker_sum[iter]
		end
	end
	# construct scope covariates
	target_liner_share = (target_liner_sum-buyer_X_scale[1]-deviater_X_scale[1])/
	                     (target_total_sum-buyer_total_sum-deviater_total_sum)
	target_special_share = (target_special_sum-buyer_X_scale[2]-deviater_X_scale[2])/
	                     (target_total_sum-buyer_total_sum-deviater_total_sum)
	target_tramper_share = (target_tramper_sum-buyer_X_scale[3]-deviater_X_scale[3])/
	                     (target_total_sum-buyer_total_sum-deviater_total_sum)
	target_tanker_share = (target_tanker_sum-buyer_X_scale[4]-deviater_X_scale[4])/
	                     (target_total_sum-buyer_total_sum-deviater_total_sum)
	# modification dropping deviator firm kk
	target_liner_sum = target_liner_sum - deviater_X_scale[1]
	target_special_sum = target_special_sum - deviater_X_scale[2]
	target_tramper_sum = target_tramper_sum - deviater_X_scale[3]
	target_tanker_sum = target_tanker_sum - deviater_X_scale[4]
	res = vcat(target_liner_sum, target_special_sum, target_tramper_sum, target_tanker_sum,
	       target_liner_share, target_special_share, target_tramper_share, target_tanker_share)
	return res
end


function extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(ii, kk, data;info_sum = temp_info_sum)
	#global liner_sum, special_sum, tanker_sum, tramper_sum
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
	group_names = data[ii,4]
	# pick individual own covariates
	buyer_X_scale = Vector(data[ii,5:8])
	buyer_total_sum = data.total[ii]
	target_total_sum = data.group_total_tonnage[ii]
	# pick unmatched kk covariates
	unmatched_X_scale = Vector(data[kk,5:8])
	unmatched_total_sum = data.total[kk]
	# construct scale covariates
	target_liner_sum = zeros(0)
	target_special_sum = zeros(0)
	target_tramper_sum = zeros(0)
	target_tanker_sum = zeros(0)
	for iter = 1:size(liner_sum,1)
		#liner_sum.group[iter] == group_names
		if liner_sum.group[iter] == group_names
			target_liner_sum = liner_sum.liner_sum[iter]
		end
		if special_sum.group[iter] == group_names
			target_special_sum = special_sum.special_sum[iter]
		end
		if tramper_sum.group[iter] == group_names
			target_tramper_sum = tramper_sum.tramper_sum[iter]
		end
		if tanker_sum.group[iter] == group_names
			target_tanker_sum = tanker_sum.tanker_sum[iter]
		end
	end
	# construct scope covariates
	target_liner_share = (target_liner_sum-buyer_X_scale[1]+unmatched_X_scale[1])/
	                     (target_total_sum-buyer_total_sum+unmatched_total_sum)
	target_special_share = (target_special_sum-buyer_X_scale[2]+unmatched_X_scale[2])/
	                     (target_total_sum-buyer_total_sum+unmatched_total_sum)
	target_tramper_share = (target_tramper_sum-buyer_X_scale[3]+unmatched_X_scale[3])/
	                     (target_total_sum-buyer_total_sum+unmatched_total_sum)
	target_tanker_share = (target_tanker_sum-buyer_X_scale[4]+unmatched_X_scale[4])/
	                     (target_total_sum-buyer_total_sum+unmatched_total_sum)
	# modification adding unmatched kk
	target_liner_sum = target_liner_sum + unmatched_X_scale[1]
	target_special_sum = target_special_sum + unmatched_X_scale[2]
	target_tramper_sum = target_tramper_sum + unmatched_X_scale[3]
	target_tanker_sum = target_tanker_sum + unmatched_X_scale[4]
	res = vcat(target_liner_sum, target_special_sum, target_tramper_sum, target_tanker_sum,
	       target_liner_share, target_special_share, target_tramper_share, target_tanker_share)
end

function gen_merger_cost_for_buyer_ii(ii, buyer1_X, target1_X, data)
	subdata = @linq data |>
		where(:group .== data[ii,4])
	# construct swapped matches
	num_of_firms_in_coalition = size(subdata, 1)
	total_tonnage = sum(buyer1_X[1:4]) + sum(target1_X[1:4])
	merger_cost = num_of_firms_in_coalition/log(total_tonnage*100 + 1) # 1million ton = 1 translated into 1ton=1
	#merger_cost = num_of_firms_in_coalition/total_tonnage
	return merger_cost
end

function gen_subsidy_indicator_for_buyer_ii(buyer1_X, target1_X,
	                                        subsidy_type;
	                                        subsidy_threshold = 1,
											subsidy_amount = 1)
	total_tonnage = sum(buyer1_X[1:4]) + sum(target1_X[1:4])
	subsidy_indicator = total_tonnage > subsidy_threshold
	if subsidy_type == "shared"
	    subsidy_effect = subsidy_amount*subsidy_indicator/total_tonnage
	elseif subsidy_type == "to_buyer"
		subsidy_effect = subsidy_amount*subsidy_indicator
	end
	return subsidy_effect
end

function gen_utility_est(ii, buyer1_X::Vector, target1_X::Vector, data,
	                     beta, gamma, delta, subsidy_type)
		#interaction_X_beta = beta.*log.(buyer1_X.*target1_X.+1)
		buyer1_X_total = sum(buyer1_X[1:4])
		target1_X_total = sum(target1_X[1:4])
		#interaction_X_beta = 1.0*buyer1_X_total.*target1_X_total .+ beta.*buyer1_X.*target1_X
		interaction_X_beta = vcat(1.0*buyer1_X_total.*target1_X_total, beta.*buyer1_X.*target1_X)
		#interaction_X_beta = beta.*buyer1_X.*target1_X
		merger_cost = gen_merger_cost_for_buyer_ii(ii, buyer1_X, target1_X, data)
		subsidy = gen_subsidy_indicator_for_buyer_ii(buyer1_X, target1_X,
		                                             subsidy_type,
		                                             subsidy_amount = 1)
		utility = sum(interaction_X_beta) - gamma*merger_cost + delta*subsidy
	return utility
end


function gen_utility_est_without_subsidy(ii, buyer1_X::Vector,
	                     target1_X::Vector, data,
	                     beta, gamma, delta)
		#interaction_X_beta = beta.*log.(buyer1_X.*target1_X.+1)
		buyer1_X_total = sum(buyer1_X[1:4])
		target1_X_total = sum(target1_X[1:4])
		interaction_X_beta = 1.0*buyer1_X_total.*target1_X_total .+ beta.*buyer1_X.*target1_X
		#interaction_X_beta = beta.*buyer1_X.*target1_X
		merger_cost = gen_merger_cost_for_buyer_ii(ii, buyer1_X, target1_X, data)
		utility = sum(interaction_X_beta) - gamma*merger_cost + 0
	return utility
end

function gen_unmatched_utility_est(buyer1_X::Vector, beta)
		#interaction_X_beta = beta.*log.(buyer1_X.+1)
		buyer1_X_total = sum(buyer1_X[1:4])
		#interaction_X_beta = 1.0*buyer1_X_total .+ beta.*buyer1_X
		#interaction_X_beta = 1.0*buyer1_X_total.*buyer1_X_total .+ beta.*buyer1_X.*buyer1_X
		interaction_X_beta = 1.0*buyer1_X_total.*0.01 .+ beta.*buyer1_X.*0.01 # specify unmatched = match with 0.01
		#interaction_X_beta = beta.*buyer1_X
		utility = sum(interaction_X_beta)
	return utility
end



#-------------------------------------#
# functions for constructing CR
# for reference
#-------------------------------------#

function inner_construct_critical_value_CR_two_variables(
	                                 subsampled_id_list,
	                                 data,
									 theta,
									 subsidy_type;
									 calibrated_delta_list = [1],
									 boot_num = 200,
									 size_of_subsample = 30,
									 info_sum = temp_info_sum)
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
	data_main_firm = @linq data |>
		where(:type .== "(1) main")
	data_not_main_firm = @linq data |>
			where(:type .!= "(1) main")
	main_firm_id = data_main_firm.firm_id
	fullsample_id = data.firm_id
	subsampled_id_list_all = zeros(Int64,
	                               boot_num,
	                               size(data_main_firm)[1] + size_of_subsample)
	for iter = 1:boot_num
		Random.seed!(iter)
	    picked_not_main_firm_id = StatsBase.sample(data_not_main_firm.firm_id,
										   size_of_subsample,
										   replace = false)
	    subsampled_id_list = vcat(main_firm_id, picked_not_main_firm_id)
		subsampled_id_list_all[iter,:] = subsampled_id_list
	end
	# construct a critical value
	single_variable_index = [1:1:length(variable_list);][variable_list .== variable][1]
	for kk = 1:length(calibrated_delta_list)
		calibrated_delta_kk = calibrated_delta_list[kk]
		function score_bthis_only_target_x(subsampled_id_list, data,
										   theta, subsidy_type;
										   calibrated_delta = calibrated_delta_kk,
										   info_sum = info_sum)
			#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
			target_theta = vcat(zeros(single_variable_index-1),
								theta[1],
								zeros(8-single_variable_index),
								theta[2],
								calibrated_delta) # first parameter must be normalized to 1
			score_res, total_num_ineq = score_b_est_data(subsampled_id_list,
			                                             data,
			                   							 target_theta,
														 subsidy_type)
			res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
			#println("score:$(res) \n" )
			return res, score_res, total_num_ineq
		end
		boot_Q = zeros(boot_num)
		# theta = [-100,10]
		@time for iter = 1:boot_num
		    boot_Q[iter] = score_bthis_only_target_x(subsampled_id_list_all[iter,:],
		                          data,
								  theta,
								  subsidy_type)[2]
		end
		full_sample_Q = score_bthis_only_target_x(fullsample_id,
							  data,
							  theta,
							  subsidy_type)[2]
	end
	return boot_Q, full_sample_Q
end

function outer_construct_critical_value_CR_two_variables(
	                                 subsampled_id_list,
	                                 data,
									 theta,
									 subsidy_type;
									 a_b = 42,
									 alpha_level = 0.05,
									 calibrated_delta_list,
									 boot_num = 200,
									 size_of_subsample = 30,
									 info_sum = temp_info_sum)
	@time boot_Q, full_sample_Q = inner_construct_critical_value_CR_two_variables(
		                                 subsampled_id_list,
		                                 data,
										 theta,
										 subsidy_type;
										 calibrated_delta_list = calibrated_delta_list,
										 boot_num = boot_num,
										 size_of_subsample = size_of_subsample,
										 info_sum = temp_info_sum)
	# find minimizer x as a critical value
	min_x = quantile(a_b.*(boot_Q), 0.93)
	x_iter = min_x
	boot_res = 0
	while boot_res <= 1 - alpha_level
		boot_res = mean(a_b.*(boot_Q) .<= x_iter)
		x_iter += 1
	end
	#@show x_iter
	return x_iter, full_sample_Q
end
# @time critical_d, full_sample_Q = outer_construct_critical_value_CR_two_variables(
# 									 subsampled_id_list,
# 									 data,
# 									 [41.1, 2.64],
# 									 subsidy_type;
# 									 a_b = 42,
# 									 alpha_level = 0.05,
# 									 calibrated_delta_list,
# 									 boot_num = 200,
# 									 size_of_subsample = 30,
# 									 info_sum = temp_info_sum)
# full_sample_Q * a_n < critical_d
# critical_d./full_sample_Q

end # end module
