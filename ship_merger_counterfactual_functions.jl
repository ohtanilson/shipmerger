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
			#interaction_X_beta_unmatched = 1.0*buyer1_X_total.*0.01 .+ beta.*buyer1_X.*0.01 # specify unmatched = match with 0.01
			interaction_X_beta_unmatched = 1.0*buyer1_X_total.*buyer1_X_total .+ beta.*buyer1_X.*buyer1_X
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

function solve_equilibrium_counterfactual(m::shipmerger_struct, utility)
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

function gen_merger_composition(m::shipmerger_struct, matches)
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

#
# function iterate_simulation_counterfactual(model_id, theta_hat_all_models,
# 	                        data,
# 	                        threshold_tonnage_list,
# 							subsidy_amount_list;
# 	                        iter_end = 1,
# 							temp_subsidy_type = "to_buyer")
# 	theta_hat = theta_hat_all_models[:,model_id]
# 	for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
# 		number_of_groups = zeros(iter_end)
# 		number_of_unmatched = zeros(iter_end)
# 		matches_list = zeros(m.N, m.PN)
# 		threshold_tonnage_iter = threshold_tonnage_list[nn]
# 		subsidy_amount_iter = subsidy_amount_list[mm]
# 	    @time for kk = 1:iter_end
# 		    utility = gen_utility_matrix_counterfactual(theta_hat, data,
# 									randomseed = kk,
# 									threshold_tonnage = threshold_tonnage_iter, # 1 million
# 									subsidy_amount = subsidy_amount_iter,
# 									subsidy_type = temp_subsidy_type)
# 		    matches = solve_equilibrium_counterfactual(m, utility)
# 			@show theta_hat
# 			@show kk, threshold_tonnage_iter, subsidy_amount_iter
# 		    @show number_of_groups[kk], number_of_unmatched[kk] = gen_merger_composition(m, matches)
# 			open("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
# 				DelimitedFiles.writedlm(io, matches,",")
# 			end
# 			println("\nthreshold_tonnage = $(threshold_tonnage_iter) and subsidy_amount = $(subsidy_amount_iter) iter = $(kk)\n")
# 			println("num of groups = $(number_of_groups[kk]) and num of unmatched = $(number_of_unmatched[kk])\n")
# 	    end
# 		open("julia_merger_result/counterfactual_number_of_groups_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
# 			DelimitedFiles.writedlm(io, number_of_groups,",")
# 		end
# 		open("julia_merger_result/counterfactual_number_of_unmatched_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
# 			DelimitedFiles.writedlm(io, number_of_unmatched,",")
# 		end
# 	end
# 	#return number_of_groups, number_of_unmatched
# end
