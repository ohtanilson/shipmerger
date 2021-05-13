using Plots
using DelimitedFiles
using Combinatorics
using CSV
using DataFrames
using DataFramesMeta
using StatsBase
using Random
using BlackBoxOptim
using LaTeXTabulars
using LaTeXStrings
using Statistics
include(joinpath(dirname(@__FILE__),"functions.jl"))
using .functions
using JuMP, Gurobi
using Distributions
using JLD
#------------------------------#
#  Counterfactual
#------------------------------#
data = CSV.read("data_for_maximum_rank_estimation.csv",DataFrame)
rename!(data, :Column1 => :firm_id)
data_for_counterfactual = @linq data |>
	where(:type .== "(1) main")
# assign theta_hat
model_specification = "two_variables_main_firms_only"
if model_specification == "column 3"
	# full model
    myests_point_full_X_only = readdlm("julia_merger_result/myests_subsample_size_106_full_X_only.txt",',',Float64)
	num_correct_ineq_full_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_full_X_only.txt",',',Float64)
	num_total_ineq_full_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_full_X_only.txt",',',Float64)
	accuracy_full_X_only = vec(num_correct_ineq_full_X_only./num_total_ineq_full_X_only)
	final_ests_point_full_X_only = round.(myests_point_full_X_only[findmax(accuracy_full_X_only)[2],:],digits=2)
	size_of_subsample_temp = 60
	myests_CI_full_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_full_X_only.txt",',',Float64)
	theta_hat_all_models = round.(hcat(vcat(final_ests_point_full_X_only[1:8],
	                final_ests_point_full_X_only[9],
					final_ests_point_full_X_only[10]),
	                # the cheapest scenario
	                vcat(final_ests_point_full_X_only[1:8],
					Statistics.quantile(myests_CI_full_X_only[:,9], [0.975])[1],
					final_ests_point_full_X_only[10]),
					# the most expensive scenario
					vcat(final_ests_point_full_X_only[1:8],
					Statistics.quantile(myests_CI_full_X_only[:,9], [0.025])[1],
					final_ests_point_full_X_only[10])
					),digits = 2)
elseif model_specification == "column 2"
	myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_106_scope_X_only.txt",',',Float64)
	num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_scope_X_only.txt",',',Float64)
	num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_scope_X_only.txt",',',Float64)
	accuracy_scope_X_only = vec(num_correct_ineq_scope_X_only./num_total_ineq_scope_X_only)
	final_ests_point_scope_X_only = round.(myests_point_scope_X_only[findmax(accuracy_scope_X_only)[2],:],digits=2)
	size_of_subsample_temp = 60
	myests_CI_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_scope_X_only.txt",',',Float64)

	theta_hat_all_models = round.(hcat(vcat(zeros(4),
	                final_ests_point_scope_X_only[1:4],
					final_ests_point_scope_X_only[5],
					final_ests_point_scope_X_only[6]),
					# the cheapest scenario
					vcat(zeros(4),
					final_ests_point_scope_X_only[1:4],
					#Statistics.quantile.(myests_CI_scope_X_only[:,1:4], [0.975])[1]
					Statistics.quantile(myests_CI_scope_X_only[:,5], [0.975])[1],
					final_ests_point_scope_X_only[6]),
					# the most expensive scenario
					vcat(zeros(4),final_ests_point_scope_X_only[1:4],
					Statistics.quantile(myests_CI_scope_X_only[:,5], [0.025])[1],
					final_ests_point_scope_X_only[6])
					),digits = 2)
elseif model_specification == "column_8_two_variables"
	# check behavior
	temp_subsidy_type = "shared"
	file_name_variable = "x8"
	size_of_subsample_temp = 30
	calibrated_delta = 1
	myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	accuracy_scope_X_only = vec(num_correct_ineq_scope_X_only./num_total_ineq_scope_X_only)
	final_ests_point_scope_X_only = round.(myests_point_scope_X_only[findmax(accuracy_scope_X_only)[2],:],digits=2)
	myests_CI_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	theta_hat_all_models = round.(hcat(vcat(zeros(7),
	                final_ests_point_scope_X_only[1],
					final_ests_point_scope_X_only[2],
					calibrated_delta),
					# the cheapest scenario
					vcat(zeros(7),
					#Statistics.quantile.(myests_CI_scope_X_only[:,1:4], [0.975])[1]
					final_ests_point_scope_X_only[1],
					Statistics.quantile(myests_CI_scope_X_only[:,2], [0.975])[1],
					calibrated_delta),
					# the most expensive scenario
					vcat(zeros(7),
					final_ests_point_scope_X_only[1],
					Statistics.quantile(myests_CI_scope_X_only[:,2], [0.025])[1],
					calibrated_delta)
					),digits = 2)
elseif model_specification == "two_variables_main_firms_only"
	# check behavior
	temp_subsidy_type = "shared"
	file_num = 1
	file_name_variable_list = ["x1","x2","x3","x4","x5","x6","x7","x8"]
	file_name_variable = file_name_variable_list[file_num]
	size_of_subsample_temp = 0
	calibrated_delta = 1000
	myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	accuracy_scope_X_only = vec(num_correct_ineq_scope_X_only./num_total_ineq_scope_X_only)
	final_ests_point_scope_X_only = round.(myests_point_scope_X_only[findmax(accuracy_scope_X_only)[2],:],digits=2)
	#myests_CI_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
	temp_myests_point_all = myests_point_scope_X_only[num_correct_ineq_scope_X_only[:,1].==Int(findmax(num_correct_ineq_scope_X_only[:,1])[1]),:]
	CI_all_table = round.(hcat(
		Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
		Statistics.quantile(temp_myests_point_all[:,2], [0.0,1.0])
		),digits=1)
	theta_hat_all_models = round.(hcat(vcat(zeros(file_num-1),
	                Statistics.quantile(temp_myests_point_all[:,1], [0.5]), # x2
					zeros(7-(file_num-1)),
					Statistics.quantile(temp_myests_point_all[:,2], [0.5]),
					calibrated_delta),
					# the cheapest scenario
					vcat(zeros(file_num-1),
					#Statistics.quantile.(myests_CI_scope_X_only[:,1:4], [0.975])[1]
					CI_all_table[2,1], # x3
					zeros(7-(file_num-1)),
					CI_all_table[1,2],
					calibrated_delta),
					# the most expensive scenario
					vcat(zeros(file_num-1),
					CI_all_table[1,1], # x3
					zeros(7-(file_num-1)),
					CI_all_table[2,2],
					calibrated_delta)
					),digits = 2)
end

# function extract_buyer_covariates_if_ii_is_buyer(ii, data)
# 	# pick individual own covariates
# 	buyer_X_scale = Vector(data[ii,5:8])
# 	buyer_X_scope = Vector(data[ii,11:14])
# 	res = vcat(buyer_X_scale,buyer_X_scope)
# 	return res
# end
#
# function gen_merger_cost_for_buyer_ii(ii, buyer1_X, target1_X, data)
# 	subdata = @linq data |>
# 		where(:group .== data[ii,4])
# 	# construct swapped matches
# 	num_of_firms_in_coalition = size(subdata, 1)
# 	total_tonnage = sum(buyer1_X[1:4]) + sum(target1_X[1:4])
# 	merger_cost = num_of_firms_in_coalition/log(total_tonnage*100 + 1) # 1million ton = 1 translated into 1ton=1
# 	#merger_cost = num_of_firms_in_coalition/total_tonnage
# 	return merger_cost
# end
#
# function gen_subsidy_indicator_for_buyer_ii(buyer1_X, target1_X,
# 	                                        subsidy_type;
# 	                                        subsidy_threshold = 1,
# 											subsidy_amount = 1)
# 	total_tonnage = sum(buyer1_X[1:4]) + sum(target1_X[1:4])
# 	subsidy_indicator = total_tonnage > subsidy_threshold
# 	if subsidy_type == "shared"
# 	    subsidy_effect = subsidy_amount*subsidy_indicator/total_tonnage
# 	elseif subsidy_type == "to_buyer"
# 		subsidy_effect = subsidy_amount*subsidy_indicator
# 	end
# 	return subsidy_effect
# end

@time m = carrier_mod(randomseed = 8,
					 N = 12,
					 ϵ_dist=Distributions.Normal(0.0,5.0),
					 β_0 = 3.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 5.0,# coefficient of additional merger cost
					 ton_dim= 2)
function gen_utility_matrix_counterfactual(theta_hat::Vector,
	                        data::DataFrame;
							randomseed = temp_randomseed,
	                        threshold_tonnage = 1, # 1 million
							subsidy_amount = 1,
							subsidy_type = "shared")
	beta = theta_hat[1:8]
	gamma = theta_hat[9]
	delta = theta_hat[10]
	@time m = carrier_mod(randomseed = randomseed,
	                      ϵ_dist=Distributions.Normal(0.0,5.0),
						  N = 12)
	m.ϵ_mat = m.ϵ_mat
	data_for_counterfactual = @linq data |>
		where(:type .== "(1) main")
	subsampled_id_list = data_for_counterfactual.firm_id
	all_possible_pairs = [Combinatorics.combinations(subsampled_id_list,2)...]
	index_list = Array{Int64,2}(undef, length(all_possible_pairs), 2)
	for i in 1:length(all_possible_pairs)
		index_list[i,1] = all_possible_pairs[i][1]
		index_list[i,2] = all_possible_pairs[i][2]
	end
	# construct utility matrix
	utility = zeros(m.N, m.PN)
	@inbounds for i = 1:m.N, j = 1:m.PN
		# pick up buyer covariates
		buyer1_X = extract_buyer_covariates_if_ii_is_buyer(subsampled_id_list[i],
		                                                   data)
		# pick up target covariates
		k_index = zeros(Int64,0)
		for k = 1:m.N # extract member of firm 1's coalition
			if m.Bundle[j][k]== '1'
				k_index = vcat(k_index, k)
			end
		end
		target1_X_each = zeros(length(k_index), 8)
		for kk = 1:length(k_index)
			firm_id_kk = k_index[kk]
			target1_X_each[kk,:] = extract_buyer_covariates_if_ii_is_buyer(subsampled_id_list[firm_id_kk], data)
		end
		target1_X_scale = sum(target1_X_each[:,1:4], dims = 1)
		target1_X_total = sum(target1_X_scale)
		target1_X_scope = target1_X_scale./target1_X_total
		target1_X = vec(hcat(target1_X_scale, target1_X_scope))
		# utility
		buyer1_X_total = sum(buyer1_X[1:4])
		target1_X_total = sum(target1_X[1:4])
		interaction_X_beta = vcat(1.0*buyer1_X_total.*target1_X_total, beta.*buyer1_X.*target1_X)
		#interaction_X_beta = beta.*buyer1_X.*target1_X
		merger_cost = gen_merger_cost_for_buyer_ii(subsampled_id_list[i], buyer1_X, target1_X, data)
		subsidy = gen_subsidy_indicator_for_buyer_ii(buyer1_X, target1_X,
		                                             subsidy_type,
		                                             subsidy_threshold = threshold_tonnage,
													 subsidy_amount = subsidy_amount)
		payoff_obs_match1 = sum(interaction_X_beta) - gamma*merger_cost + delta*subsidy
	    if IS_Sell(m,i,j)
		    # match with himself = (unmatched)
			buyer1_X_total = sum(buyer1_X[1:4])
			interaction_X_beta_unmatched = 1.0*buyer1_X_total.*0.01 .+ beta.*buyer1_X.*0.01 # specify unmatched = match with 0.01
		    utility[i,j] = sum(interaction_X_beta_unmatched) + m.ϵ_mat[i,j]
		else
			utility[i,j] = payoff_obs_match1 + m.ϵ_mat[i,j] # ad hoc scaling
	    end
	end
	# assigning negative infinity to unreal allocation. (ex.) (1, (101)) and (1,(1,1,1))
	@inbounds for i = 1:m.N, j = 1:m.PN
            if IS_outdomain(m.N,i,j)      # hard code an continuum agent model to reach an integer equilibrium
                utility[i,j] = -99999999
            elseif IS_Sell(m,i,j)
                utility[i,1] = copy(utility[i,j]) # unmatched payoff
				utility[i,j] = m.ϵ_mat[i,j]  # selling payoff
            else
				utility[i,j] = utility[i,j] # keep original payoff
			end
    end
	return utility
end

function solve_equilibrium_counterfactual(m::carrier_struct_v1, utility)
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
	function buyer_expression_ipopt(m::carrier_struct_v1, t::Integer)
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
				double_cc_expression_1 = x[m,n] # sell side
				double_cc_expression_2 = buyer_expression_ipopt(mm,m) # buy side
				@constraint(model,double_cc_expression_1 == double_cc_expression_2)
			else
				nothing
			end
		end
		@show m
	end
	println( "# feasibility constraints(measure const):\n")
	@constraint(model, feas_rho[i=1:N], sum(x[i,j] for j in 1:PN)== eta)

    JuMP.@objective(model, Max, sum(x[i,j]*utility[i,j] for i in 1:N, j in 1:PN))
    println("Time for optimizing model:")
    @time JuMP.optimize!(model)
    # show results
    objv = JuMP.objective_value(model)
    println("objvalue　= ", objv)
    matches = JuMP.value.(x)
	matches = round.(matches, digits = 1)
	res = copy(matches) .+ rand(m.N, m.PN).*1e-8 # add small pertubations
	for i = 1:m.N
		# to avoid real-valued solution, add randomness and choose maximum
		res[i,findmax(res[i,:])[2]] = 1
	end
	res = res.==1
    return res
end

function gen_merger_composition(m::carrier_struct_v1, matches)
	if sum(matches.==1) != m.N
		# the equilibrium is not integer
        return 0, 0
	else
		obs_bundle = zeros(Int64, m.N)
		# pick up observed merger bundle for each buyer
		for i = 1:m.N
			obs_bundle[i] = Int([1:1:m.PN;][matches[i,:] .== 1][1])
		end
		# count unique group buyer
		buyer_bundle = zeros(0)
		for i = 1:length(obs_bundle)
			is_sell_vector = IS_Sell(m, i, obs_bundle[i])
			if is_sell_vector == 0
				buyer_bundle = vcat(buyer_bundle, obs_bundle[i])
			end
		end
		number_of_groups = sum(unique(buyer_bundle).!=1)
		number_of_unmatched = sum(buyer_bundle.==1)
		return number_of_groups, number_of_unmatched
	end
end
#m.Bundle[Int.(buyer_bundle)]

function iterate_simulation_counterfactual(model_id, theta_hat_all_models,
	                        data,
	                        threshold_tonnage_list,
							subsidy_amount_list;
	                        iter_end = 1)
	theta_hat = theta_hat_all_models[:,model_id]
	for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
		number_of_groups = zeros(iter_end)
		number_of_unmatched = zeros(iter_end)
		matches_list = zeros(m.N, m.PN)
		threshold_tonnage_iter = threshold_tonnage_list[nn]
		subsidy_amount_iter = subsidy_amount_list[mm]
	    @time for kk = 1:iter_end
		    utility = gen_utility_matrix_counterfactual(theta_hat, data,
									randomseed = kk,
									threshold_tonnage = threshold_tonnage_iter, # 1 million
									subsidy_amount = subsidy_amount_iter,
									subsidy_type = "shared")
		    matches = solve_equilibrium_counterfactual(m, utility)
			@show theta_hat
			@show kk, threshold_tonnage_iter, subsidy_amount_iter
		    @show number_of_groups[kk], number_of_unmatched[kk] = gen_merger_composition(m, matches)
			open("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
				DelimitedFiles.writedlm(io, matches,",")
			end
			println("\nthreshold_tonnage = $(threshold_tonnage_iter) and subsidy_amount = $(subsidy_amount_iter) iter = $(kk)\n")
			println("num of groups = $(number_of_groups[kk]) and num of unmatched = $(number_of_unmatched[kk])\n")
	    end
		open("julia_merger_result/counterfactual_number_of_groups_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
			DelimitedFiles.writedlm(io, number_of_groups,",")
		end
		open("julia_merger_result/counterfactual_number_of_unmatched_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
			DelimitedFiles.writedlm(io, number_of_unmatched,",")
		end
	end
	#return number_of_groups, number_of_unmatched
end
maximum(data.total)
# test a single computation
temp_model_id = 3#2
theta_hat = theta_hat_all_models[:,temp_model_id]
#theta_hat[9] = theta_hat[9]*0.01
#theta_hat[1:4] .= 0 # drop scale covariates
positive_constant = 180
while number_of_groups_temp < 6
	# find positive constant term for correcting unmatched and matched
	@time utility = gen_utility_matrix_counterfactual(theta_hat,
	                        data,
							randomseed = 1,
							threshold_tonnage = 1, # 1 million
							subsidy_amount = 1,#0,
							subsidy_type = "shared")
	modified_utility = copy(utility)
	for i = 1:m.N
		for j = 2:m.PN
			if IS_Sell(m,i,j)!=1
		        modified_utility[i,j] = utility[i,j] - positive_constant
			end
		end
	end
	matches = solve_equilibrium_counterfactual(m, modified_utility)
	#matches = solve_equilibrium_counterfactual(m, utility)
	for ii = 1:m.N
		@show ii
		@show findmax(matches[ii,:].>0)
		@show modified_utility[ii,matches[ii,:]]
	    @show m.Bundle[matches[ii,:].>0]
	end
    @show number_of_groups_temp, number_of_unmatched_temp = gen_merger_composition(m, matches)
	@show positive_constant += 1
	# positive_constant = 240
end

threshold_tonnage_list = [1, 2, 3, 4, 5, 7.5]# 6dim
subsidy_amount_list = [0, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 5] # 6 dim (->9dim)
model_list = [1, 2, 3]
iter_end = 1
Threads.nthreads()
Threads.threadid()
JULIA_NUM_THREADS=8
#positive_constant = 140 #135
@time Threads.@threads for model_id_iter in model_list
	theta_hat = theta_hat_all_models[:,model_id_iter]
    for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
		number_of_groups = zeros(iter_end)
		number_of_unmatched = zeros(iter_end)
		matches_list = zeros(m.N, m.PN)
		threshold_tonnage_iter = threshold_tonnage_list[nn]
		subsidy_amount_iter = subsidy_amount_list[mm]
	    @time for kk = 1:iter_end
		    utility = gen_utility_matrix_counterfactual(theta_hat, data,
									randomseed = kk,
									threshold_tonnage = threshold_tonnage_iter, # 1 million
									subsidy_amount = subsidy_amount_iter,
									subsidy_type = "shared")
			modified_utility = copy(utility)
			# correction constant term for payoffs of being matched
			if model_id_iter == 2
				# cheapest scenario
				positive_constant = 180
			else
				# the most expensive and middle scenario
				positive_constant = 140
			end
			for i = 1:m.N
				for j = 2:m.PN
					if IS_Sell(m,i,j)!=1
				        modified_utility[i,j] = utility[i,j] - positive_constant
					end
				end
			end
		    matches = solve_equilibrium_counterfactual(m, modified_utility)
			for ii = 1:m.N
				@show ii
				@show findmax(matches[ii,:].>0)
				@show modified_utility[ii,matches[ii,:]]
			    @show m.Bundle[matches[ii,:].>0]
			end
			@show theta_hat
			@show kk, threshold_tonnage_iter, subsidy_amount_iter
		    @show number_of_groups[kk], number_of_unmatched[kk] = gen_merger_composition(m, matches)
			open("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
				DelimitedFiles.writedlm(io, convert(Array{Float64,2},matches),",")
			end
			println("\nthreshold_tonnage = $(threshold_tonnage_iter) and subsidy_amount = $(subsidy_amount_iter) iter = $(kk)\n")
			println("num of groups = $(number_of_groups[kk]) and num of unmatched = $(number_of_unmatched[kk])\n")
	    end
		open("julia_merger_result/counterfactual_number_of_groups_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
			DelimitedFiles.writedlm(io, number_of_groups,",")
		end
		open("julia_merger_result/counterfactual_number_of_unmatched_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
			DelimitedFiles.writedlm(io, number_of_unmatched,",")
		end
	end
	# 6dim*6dim*3scenarios
	#36585.071584 seconds (89.23 G allocations: 26.286 TiB, 32.74% gc time, 0.00% compilation time)
	# 6dim*6dim*3scenarios on 210505
	# 42002.569748 seconds (89.22 G allocations: 26.286 TiB, 32.46% gc time, 0.00% compilation time)
	# 12601.544318 seconds (89.22 G allocations: 26.286 TiB, 27.80% gc time, 0.00% compilation time)
	# 6dim*9dim*3scenarios on 210506
	#27590.268983 seconds (133.92 G allocations: 39.436 TiB, 46.92% gc time, 0.00% compilation time
end
kk = 1
model_id_iter = 1
threshold_tonnage_iter = 1.0
subsidy_amount_iter = 0.0
temp_d = readdlm("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
temp_iter_kk = 1
using LaTeXTabulars
using LaTeXStrings
num_group_list = zeros(length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list))
num_unmatched_list = zeros(length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list))
matches_list = zeros(length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list),
					  m.N,
					  m.PN)
@time for model_id in model_list
	for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
		threshold_tonnage_iter = threshold_tonnage_list[nn]
		subsidy_amount_iter = subsidy_amount_list[mm]
		counterfactual_number_of_groups = readdlm("julia_merger_result/counterfactual_number_of_groups_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
		counterfactual_number_of_unmatched = readdlm("julia_merger_result/counterfactual_number_of_unmatched_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
		num_group_list[model_id, nn, mm] = median(counterfactual_number_of_groups)
		num_unmatched_list[model_id, nn, mm] = median(counterfactual_number_of_unmatched)
		println("model id = $(model_id) and threshold_tonnage = $(threshold_tonnage_iter) and subsidy_amount = $(subsidy_amount_iter)")
		println("num_group=",num_group_list[model_id, nn, mm],
		        ", unmatched=",num_unmatched_list[model_id, nn, mm])
	end
	for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
		threshold_tonnage_iter = threshold_tonnage_list[nn]
		subsidy_amount_iter = subsidy_amount_list[mm]
		counterfactual_matches = zeros(m.N, m.PN, iter_end)
		for kk = 1:iter_end
	        counterfactual_matches[:,:,kk] = readdlm("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
		end
		matches_list[model_id, nn, mm, :, :] = mean(counterfactual_matches, dims = 3)
	end
end


# LaTeXTabulars.latex_tabular("julia_merger_table/counterfactual_subsidy_threshold.tex",
#               Tabular("@{\\extracolsep{5pt}}lcccc"),
#               [Rule(:top),
#                ["","","Point Estimate Scenario", "The most expensive scenario", "The cheapest expenditure scenario"],
# 			   ["merger cost \$\\gamma\$","",
# 			   "\$\\gamma=-$(theta_hat_all_models[9,1])\$",
# 			   "\$\\gamma=-$(theta_hat_all_models[9,2])\$",
# 			   "\$\\gamma=-$(theta_hat_all_models[9,3])\$"],
# 			   ["","",
# 			    "Avg \$\\sharp\$ of groups (\$\\sharp\$ of unmatched)",
# 			    "Avg \$\\sharp\$ of groups (\$\\sharp\$ of unmatched)",
# 				"Avg \$\\sharp\$ of groups (\$\\sharp\$ of unmatched)"],
# 			   Rule(:mid),
# 			   ["subsidy threshold (DW)", "", "", "", ""],
# 			   #beta_0
# 			   ["", "", "", "", ""],
# 			   ["$(threshold_tonnage_list[1]) million ton (benchmark, data)", "",
# 			   "$(num_group_list[1,1,2]) ($(num_unmatched_list[1,1,2]))",
# 			   "$(num_group_list[2,1,2]) ($(num_unmatched_list[2,1,2]))",
# 			   "$(num_group_list[3,1,2]) ($(num_unmatched_list[3,1,2]))"],
# 			   ["" , "" , "", "", ""],
# 			   #beta_1
# 			   ["$(threshold_tonnage_list[2]) million ton", "",
# 			   "$(num_group_list[1,2,2]) ($(num_unmatched_list[1,2,2]))",
# 			   "$(num_group_list[2,2,2]) ($(num_unmatched_list[2,2,2]))",
# 			   "$(num_group_list[3,2,2]) ($(num_unmatched_list[3,2,2]))"],
# 			   ["" , "" , "", "", ""],
# 			   ["$(threshold_tonnage_list[3]) million ton" , "" ,
# 			   "$(num_group_list[1,3,2]) ($(num_unmatched_list[1,3,2]))",
# 			   "$(num_group_list[2,3,2]) ($(num_unmatched_list[2,3,2]))",
# 			   "$(num_group_list[3,3,2]) ($(num_unmatched_list[3,3,2]))"],
# 			   ["" , "" , "", "", ""],
# 			   #beta_2
#                ["$(threshold_tonnage_list[4]) million ton", "",
# 			   "$(num_group_list[1,4,2]) ($(num_unmatched_list[1,4,2]))",
# 			   "$(num_group_list[2,4,2]) ($(num_unmatched_list[2,4,2]))",
# 			   "$(num_group_list[3,4,2]) ($(num_unmatched_list[3,4,2]))"],
# 			   ["" , "" , "", "", ""],
# 			   ["$(threshold_tonnage_list[5]) million ton" , "" ,
# 			   "$(num_group_list[1,5,2]) ($(num_unmatched_list[1,5,2]))",
# 			   "$(num_group_list[2,5,2]) ($(num_unmatched_list[2,5,2]))",
# 			   "$(num_group_list[3,5,2]) ($(num_unmatched_list[3,5,2]))"],
# 			   #beta_3
# 			   ["" , "" , "", "", ""],
# 			   ["$(threshold_tonnage_list[6]) million ton" , "" ,
# 			   "$(num_group_list[1,6,2]) ($(num_unmatched_list[1,6,2]))",
# 			   "$(num_group_list[2,6,2]) ($(num_unmatched_list[2,6,2]))",
# 			   "$(num_group_list[3,6,2]) ($(num_unmatched_list[3,6,2]))"],
# 			   # ["" , "" , "", "", ""],
# 			   # ["$(threshold_tonnage_list[7]) million ton" , "" ,
# 			   # "$(num_group_list[1,7,2]) ($(num_unmatched_list[1,7,2]))",
# 			   # "$(num_group_list[2,7,2]) ($(num_unmatched_list[2,7,2]))",
# 			   # "$(num_group_list[3,7,2]) ($(num_unmatched_list[3,7,2]))"],
#                Rule(),           # a nice \hline to make it ugly
#                Rule(:bottom)])
#
#
# LaTeXTabulars.latex_tabular("julia_merger_table/counterfactual_subsidy_amount.tex",
#               Tabular("@{\\extracolsep{5pt}}lcccc"),
#               [Rule(:top),
#                ["","","Point Estimate Scenario", "The most expensive scenario", "The cheapest expenditure scenario"],
# 			   ["merger cost \$\\gamma\$","",
# 			   "\$\\gamma=-$(theta_hat_all_models[9,1])\$",
# 			   "\$\\gamma=-$(theta_hat_all_models[9,2])\$",
# 			   "\$\\gamma=-$(theta_hat_all_models[9,3])\$"],
# 			   ["","",
# 			    "Avg \$\\sharp\$ of groups (\$\\sharp\$ of unmatched)",
# 			    "Avg \$\\sharp\$ of groups (\$\\sharp\$ of unmatched)",
# 				"Avg \$\\sharp\$ of groups (\$\\sharp\$ of unmatched)"],
# 			   Rule(:mid),
# 			   ["subsidy amount (\$M\$)", "", "", "", ""],
# 			   #beta_0
# 			   ["", "", "", "", ""],
# 			   ["$(subsidy_amount_list[1])", "",
# 			   "$(num_group_list[1,1,1]) ($(num_unmatched_list[1,1,1]))",
# 			   "$(num_group_list[2,1,1]) ($(num_unmatched_list[2,1,1]))",
# 			   "$(num_group_list[3,1,1]) ($(num_unmatched_list[3,1,1]))"],
# 			   ["" , "" , "", "", ""],
# 			   #beta_1
# 			   ["$(subsidy_amount_list[2])(benchmark, data)", "",
# 			   "$(num_group_list[1,1,2]) ($(num_unmatched_list[1,1,2]))",
# 			   "$(num_group_list[2,1,2]) ($(num_unmatched_list[2,1,2]))",
# 			   "$(num_group_list[3,1,2]) ($(num_unmatched_list[3,1,2]))"],
# 			   ["" , "" , "", "", ""],
# 			   ["$(subsidy_amount_list[3])" , "" ,
# 			   "$(num_group_list[1,1,3]) ($(num_unmatched_list[1,1,3]))",
# 			   "$(num_group_list[2,1,3]) ($(num_unmatched_list[2,1,3]))",
# 			   "$(num_group_list[3,1,3]) ($(num_unmatched_list[3,1,3]))"],
# 			   ["" , "" , "", "", ""],
# 			   #beta_2
#                ["$(subsidy_amount_list[4])", "",
# 			   "$(num_group_list[1,1,4]) ($(num_unmatched_list[1,1,4]))",
# 			   "$(num_group_list[2,1,4]) ($(num_unmatched_list[2,1,4]))",
# 			   "$(num_group_list[3,1,4]) ($(num_unmatched_list[3,1,4]))"],
# 			   ["" , "" , "", "", ""],
# 			   ["$(subsidy_amount_list[5])" , "" ,
# 			   "$(num_group_list[1,1,5]) ($(num_unmatched_list[1,1,5]))",
# 			   "$(num_group_list[2,1,5]) ($(num_unmatched_list[2,1,5]))",
# 			   "$(num_group_list[3,1,5]) ($(num_unmatched_list[3,1,5]))"],
# 			   #beta_3
# 			   ["" , "" , "", "", ""],
# 			   ["$(subsidy_amount_list[6])" , "" ,
# 			   "$(num_group_list[1,1,6]) ($(num_unmatched_list[1,1,6]))",
# 			   "$(num_group_list[2,1,6]) ($(num_unmatched_list[2,1,6]))",
# 			   "$(num_group_list[3,1,6]) ($(num_unmatched_list[3,1,6]))"],
#                Rule(),           # a nice \hline to make it ugly
#                Rule(:bottom)])

Plots.plot(title = "Number of groups (δ = $(theta_hat_all_models[10,1])), Middle scenario",
           xlabel = "Subsidy amount")
Plots.vline!([1], color = :black, style = :dash, label = "")
Plots.hline!([6], color = :black, style = :dash, label = "(Data)")
for ss = 1:length(num_group_list[1,:,1])
    Plots.plot!(subsidy_amount_list,
	           num_group_list[1,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss])",
			   markershape = :auto,
			   alpha = 0.6)
end
Plots.plot!()
savefig("julia_merger_figure/counterfactual_num_of_groups_$(temp_subsidy_type)_subsidy")

Plots.plot(title = "Number of unmatched firms (δ = $(theta_hat_all_models[10,1])), Middle scenario",
           xlabel = "Subsidy amount",
		   ylim = [0, 12.5])
Plots.vline!([1], color = :black, style = :dash, label = "")
Plots.hline!([0], color = :black, style = :dash, label = "(Data)")
for ss = 1:length(num_unmatched_list[1,:,1])
    Plots.plot!(subsidy_amount_list,
	           num_unmatched_list[1,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss])",
			   markershape = :auto,
			   alpha = 0.6)
end
Plots.plot!()
savefig("julia_merger_figure/counterfactual_num_of_unmatched_$(temp_subsidy_type)_subsidy")

@show num_group_list[1,:,:].-num_group_list[2,:,:]
@show num_unmatched_list[1,:,:].-num_unmatched_list[2,:,:]
@show num_unmatched_list[1,:,:].-num_unmatched_list[3,:,:]
@show num_group_list[1,:,:].-num_group_list[3,:,:]

Plots.plot(title = "Number of groups (δ = $(theta_hat_all_models[10,1]))",
           xlabel = "Subsidy amount",
		   ylim = [0, 7])
Plots.vline!([1], color = :black, style = :dash, label = "")
for ss = 1:3
    Plots.plot!(subsidy_amount_list,
	           num_group_list[1,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Middle",
			   markershape = :auto,
			   alpha = 0.6)
end
for ss = 1:3
    Plots.plot!(subsidy_amount_list,
	           num_group_list[2,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Cheapest",
			   style = :dash,
			   alpha = 0.6)
end
for ss = 1:3
    Plots.plot!(subsidy_amount_list,
	           num_group_list[3,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Expensive",
			   style = :dot,
			   alpha = 0.6)
end
Plots.plot!()
savefig("julia_merger_figure/counterfactual_num_of_groups_different_scenarios_$(temp_subsidy_type)_subsidy")

Plots.plot(title = "Number of unmatched firms (δ = $(theta_hat_all_models[10,1]))",
           xlabel = "Subsidy amount",
		   ylim = [0, 12.5])
Plots.vline!([1], color = :black, style = :dash, label = "")
for ss = 1:3
    Plots.plot!(subsidy_amount_list,
	           num_unmatched_list[1,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Middle",
			   markershape = :auto,
			   alpha = 0.6)
end
for ss = 1:3
    Plots.plot!(subsidy_amount_list,
	           num_unmatched_list[2,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Cheapest",
			   style = :dash,
			   alpha = 0.6)
end
for ss = 1:3
    Plots.plot!(subsidy_amount_list,
	           num_unmatched_list[3,ss,:],
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Expensive",
			   style = :dot,
			   alpha = 0.6)
end
Plots.plot!()
savefig("julia_merger_figure/counterfactual_num_of_unmatched_different_scenarios_$(temp_subsidy_type)_subsidy")


Plots.plot(title = "Total Expenditure (amount × num of groups) (δ = $(theta_hat_all_models[10,1]))",
           ylabel = "Total Expenditure",
		   xlabel = "Subsidy amount",
		   ylim = [0, 24])
Plots.vline!([1], color = :black, style = :dash, label = "")
Plots.hline!([6], color = :black, style = :dash, label = "Expenditure (Data)")
for ss = 1:1:length(num_group_list[1,:,1])
    Plots.plot!(subsidy_amount_list,
	           num_group_list[1,ss,:].*subsidy_amount_list,
	           label = "subsidy threshold = $(threshold_tonnage_list[ss]), Middle",
			   markershape = :auto,
			   alpha = 0.6)
end
Plots.plot!()
savefig("julia_merger_figure/counterfactual_total_expenditure_$(temp_subsidy_type)_subsidy")

# For sankey diagram
model_id_iter = 1 # middle scenario
threshold_tonnage_list = [1, 2, 3, 4, 5, 7.5]# 6dim
subsidy_amount_list = [0, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 5] # 6 dim (->9dim)
model_list = [1, 2, 3]
target_matches_list = matches_list[model_id_iter,:,:,:,:]

# Case 1: actual amount with different thresholds
temp_matches = target_matches_list[:,6,:,:]
for temp_threshold = 1:length(threshold_tonnage_list)
	println("*****$(threshold_tonnage_list[temp_threshold])******************")
	#@show threshold_tonnage_list[temp_threshold]
	for ii = 1:m.N
		#@show ii
		#@show findmax(temp_matches[temp_threshold,ii,:].>0)
		#@show modified_utility[ii,matches[ii,:]]
		m.Bundle[temp_matches[temp_threshold,ii,:].>0]
		if m.Bundle[temp_matches[temp_threshold,ii,:].>0] != "000000000000"
		    println("firm $ii buys $(m.Bundle[temp_matches[temp_threshold,ii,:].>=1.0])")
		else
			println("firm $ii is unmatched")
		end
	end
end
println("
threshold_tonnage_list[temp_threshold] = 1.0
1 buys 4
5 buys 10
6 buys 3
9 buys 2
11 buys 7
12 buys 8
")

println("
threshold_tonnage_list[temp_threshold] = 2.0
1 buys 4 (8 typo)
5 buys 7,10 (9,12 typo?)
6 buys 2,3,11,
89,12 unmatched (typo?)
")

println("
threshold_tonnage_list[temp_threshold] = 3.0
1 buys 7,10 (8 typo?)
5 buys 4,9,12
6 buys 2,3,11
")

println("
threshold_tonnage_list[temp_threshold] = 4.0
1 buys 3,4,7,9,10(typo?)
6 buys 2,5,8,10,11
12 unmatched
")

println("
threshold_tonnage_list[temp_threshold] = 5.0
1 buys 2,3,4,6,7,8,9,10
5 unmatched
11 unmatched
12 unmatched
")

println("
threshold_tonnage_list[temp_threshold] = 7.5
1 buys 2,3,4,5,6,7,8,9,10,11,12
")


# Case 2: different amount with an actual threshold
temp_matches = target_matches_list[1,:,:,:]
for temp_amount = 1:length(subsidy_amount_list)
	println("******amount = $(subsidy_amount_list[temp_amount])*****************")
	subsidy_amount_list[temp_amount]
	for ii = 1:m.N
		#@show ii
		#@show findmax(temp_matches[temp_threshold,ii,:].>0)
		#@show modified_utility[ii,matches[ii,:]]
		m.Bundle[temp_matches[temp_amount,ii,:].>0]
		if m.Bundle[temp_matches[temp_amount,ii,:].>0] != "000000000000"
		    println("firm $ii buys $(m.Bundle[temp_matches[temp_amount,ii,:].>0])")
		else
			println("firm $ii is unmatched")
		end
	end
end
println("
subsidy_amount_list[temp_amount] = 0.0
all unmatched
")
println("
subsidy_amount_list[temp_amount] = 0.1
1 buys 23456789,10,12
11 unmatched
")
println("
subsidy_amount_list[temp_amount] = 0.25
1 buys 23456789,10,12
11 unmatched
")
println("
subsidy_amount_list[temp_amount] = 0.5
1 unmatched
4 unmatched
5 buys 10
6 buys 3
9 buys 2
11 buys 7
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 0.75
1 buys 4
5 buys 10
6 buys 3
9 buys 2
11 buys 7
12 buys 8
")
println("
subsidy_amount_list[temp_amount] = 1.0
1 buys 4
5 buys 10
6 buys 3
9 buys 2
11 buys 7
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 2.0
1 buys 4
5 buys 10
6 buys 7
9 buys 2
11 buys 3
12 buys 8
")

# end script
