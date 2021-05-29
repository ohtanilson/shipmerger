println("This program solves one-to-many one-sided TU matching problem")
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
temp_randomseed = 1
#-----------------------#
# comparative statics
#-----------------------#
@time m = carrier_mod(randomseed = 8,
					 N = 8,
					 β_0 = 3.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 5.0,# coefficient of additional merger cost
					 ton_dim= 2)

#------------------------------#
#  Estimation
#------------------------------#
data = CSV.read("data_for_maximum_rank_estimation.csv",DataFrame)
rename!(data, :Column1 => :firm_id)

#num_agents = size(data)[1]
names(data)
gdf = groupby(data, :group)
#=
liner_sum = combine(gdf, [:liner] => sum)
special_sum = combine(gdf, [:special] => sum)
tanker_sum = combine(gdf, [:tanker] => sum)
tramper_sum = combine(gdf, [:tramper] => sum)
=#
temp_info_sum = Dict([("liner_sum", combine(gdf, [:liner] => sum)),
                 ("special_sum", combine(gdf, [:special] => sum)),
				 ("tanker_sum", combine(gdf, [:tanker] => sum)),
				 ("tramper_sum", combine(gdf, [:tramper] => sum))])



function score_b_est_data(subsampled_id_list,
	                      data,
	                      theta::Vector{Float64},
						  subsidy_type;
						  info_sum = temp_info_sum,
						  compare_matched_and_unmatched = "no",
						  compare_with_and_without_subsidy = "yes")
	#global liner_sum, special_sum, tanker_sum, tramper_sum
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]

	beta = theta[1:8]
	gamma = theta[9] # coefficient on merger cost
	delta = theta[10] # coefficient on subsidy indicator
	# pick up an agent pair
	all_possible_pairs = [Combinatorics.combinations(subsampled_id_list,2)...]
	index_list = Array{Int64,2}(undef, length(all_possible_pairs), 2)
	for i in 1:length(all_possible_pairs)
		index_list[i,1] = all_possible_pairs[i][1]
		index_list[i,2] = all_possible_pairs[i][2]
	end
	# identify the role of firms (buyer, seller, unmatched)
	type_list_firm1 = data[index_list[:,1],3]
	type_list_firm2 = data[index_list[:,2],3]
	ineq = zeros(length(index_list[:,1]), 100)
	# pick up buyer covariates
	@inbounds for i = 1:length(index_list[:,1])
		idx = index_list[i,:] # pick buyer id pair
		#target_bundle_id = [matching_index[idx[1],2], matching_index[idx[2],2]]
		if type_list_firm1[i] == "(1) main" && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 1: both firms are buyers.")
			# pick up buyer covariates
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			buyer2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum)
			# pick up target covariates
		    target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data, info_sum = temp_info_sum)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum)
			#First, I construct matching maximum score inequality for two acquiring firms without price data:
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target2_X, data, beta, gamma, delta, subsidy_type) # swapped
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target1_X, data, beta, gamma, delta, subsidy_type) # swapped
			#ineq[i, kk, hh] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			ineq[i, 1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				# buyer 1
				payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, 2] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
				ineq[i, 3] = payoff_obs_match1 - payoff_obs_unmatch1 # with subsidy
				# buyer 2
				payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, 4] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
				ineq[i, 5] = payoff_obs_match2 - payoff_obs_unmatch2 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, 2] = payoff_obs_match1 - payoff_obs_unmatch1
				ineq[i, 3] = payoff_obs_match2 - payoff_obs_unmatch2
			end
		#elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && IS_Sell(m,idx[2],target_bundle_id[2])
		elseif type_list_firm1[i] == "(1) main" && (type_list_firm2[i] != "(1) main" && type_list_firm2[i] != "unmatched")
			#println("iter $i = Case 2: firm 1 is a buyer and firm 2 is a seller.")
			#Second, I construct inequalities from an observed coalition:
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data, info_sum = temp_info_sum)
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition of buyer's bundle
			subdata = @linq data |>
			    where(:group .== data[idx[1],4])
		    # construct swapped matches
			for kk = 1:size(subdata, 1)
				if subdata.id[kk] != data.id[idx[1]]
					target1_X_dropping_kk = extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(idx[1], kk, data, info_sum = temp_info_sum)
					payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X_dropping_kk, data, beta, gamma, delta, subsidy_type)
					buyer2_X_unmatched_deviator_kk = extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(idx[1], kk, data)
					payoff_unobs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched_deviator_kk, beta)
					ineq[i,kk] = payoff_obs_match1 - payoff_unobs_match1 - payoff_unobs_match2
				else
					ineq[i,kk] = 0 # buyer and dropped firm are the same so ill-defined
					#println("buyer and dropped firm are the same so null valued.")
				end
			end
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				# buyer 1
				payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, size(subdata, 1)+1] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
				ineq[i, size(subdata, 1)+2] = payoff_obs_match1 - payoff_obs_unmatch1 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, size(subdata, 1)+3] = payoff_obs_match1 - payoff_obs_unmatch1
			end
		#elseif IS_Sell(m,idx[1],target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
	    elseif (type_list_firm1[i] != "(1) main" && type_list_firm1[i] != "unmatched") && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 3: firm 1 is a seller and firm 2 is a buyer.")
			#Second, I construct inequalities from an observed coalition:
			buyer2_X = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum)
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition of buyer's bundle
			subdata = @linq data |>
			    where(:group .== data[idx[2],4])
			# construct swapped matches
			for kk = 1:size(subdata, 1)
				if subdata.id[kk] != data.id[idx[2]]
					target2_X_dropping_kk = extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(idx[2], kk, data, info_sum = temp_info_sum)
					payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X_dropping_kk, data, beta, gamma, delta, subsidy_type)
					buyer1_X_unmatched_deviator_kk = extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(idx[2], kk, data)
					payoff_unobs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched_deviator_kk, beta)
					ineq[i,kk] = payoff_obs_match2 - payoff_unobs_match2 - payoff_unobs_match1
				else
				    ineq[i,kk] = 0 # buyer and dropped firm are the same so ill-defined
					#println("buyer and dropped firm are the same so null valued.")
				end
			end
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				# buyer 2
				payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, size(subdata, 1)+1] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
				ineq[i, size(subdata, 1)+2] = payoff_obs_match2 - payoff_obs_unmatch2 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, size(subdata, 1)+3] = payoff_obs_match2 - payoff_obs_unmatch2
			end
		#elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
	    elseif type_list_firm1[i] == "(1) main" && type_list_firm2[i] == "unmatched"
			#println("iter $i = Case 4: firm 1 is a buyer and firm 2 is unmatched.")
			#Third, I construct inequalities from an unmatched target:
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data, info_sum = temp_info_sum)
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta) # unmatched
			# construct swapped matches
			target1_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[1], idx[2], data, info_sum = temp_info_sum)
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i, 1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				# buyer 1
				payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, 2] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
				ineq[i, 3] = payoff_obs_match1 - payoff_obs_unmatch1 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, 4] = payoff_obs_match1 - payoff_obs_unmatch1
			end
		#elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
	    elseif type_list_firm1[i] == "unmatched" && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 5: firm 1 is unmatched and firm 2 is a buyer.")
			#Third, I construct inequalities from an unmatched target:
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta) # unmatched
			buyer2_X = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum)
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition
			# construct swapped matches
			target2_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[2], idx[1], data, info_sum = temp_info_sum)
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match2
			# buyer 2
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i,2] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
				ineq[i,3] = payoff_obs_match2 - payoff_obs_unmatch2 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i,4] = payoff_obs_match2 - payoff_obs_unmatch2
			end
		#elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
	    elseif type_list_firm1[i] == "unmatched" && type_list_firm2[i] == "unmatched"
			#println("iter $i = Case 6: both picked firms are unmatched.")
			# Fourth, I construct inequalities from IR conditions without subsidy
			# Here, deviation comes from merging unmatched pairwise firm
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta) # unmatched
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta) # unmatched
			# construct swapped matches
			target1_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[1], idx[2], data, info_sum = temp_info_sum)
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X_unmatched, target1_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			target2_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[2], idx[1], data, info_sum = temp_info_sum)
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X_unmatched, target2_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i,1,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
		else
			#println("iter $i = Case 7,8,9: no firms are buyers.")
			ineq[i, 1] = 0
		end
	end
	total_num_ineq = sum(ineq .!= 0)
	res = sum(ineq.>0)
	return res, total_num_ineq
end

#res,total_num_ineq = score_b_est_data(subsampled_id_list, data, theta)
function iterate_estimation(data,
	                        score_bthis::Function,
							subsidy_type
	                        ;param_dim=9,
	                        num_its=10,
	                        size_of_subsample=60,
							num_steps_DE=10,
							temp_calibrated_delta = 1,
							info_sum = temp_info_sum)
	#global liner_sum, special_sum, tanker_sum, tramper_sum
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
	myests = zeros(num_its, param_dim)
	num_correct_ineq = zeros(num_its, 1)
	num_total_ineq = zeros(num_its, 1)
	data_main_firm = @linq data |>
		where(:type .== "(1) main")
	data_not_main_firm = @linq data |>
			where(:type .!= "(1) main")
	main_firm_id = data_main_firm.firm_id
	for iter = 1:num_its
		Random.seed!(iter)
		picked_not_main_firm_id = StatsBase.sample(data_not_main_firm.firm_id,
		                                           size_of_subsample,
												   replace = false)
		subsampled_id_list = vcat(main_firm_id, picked_not_main_firm_id)
		# Estimation
		temp_search_domain = [(-200.0,200.0),(-200.0,200.0),(-200.0,200.0),(-200.0,200.0),
		                	  (-200.0,200.0),(-200.0,200.0),(-200.0,200.0),(-200.0,200.0)]
		m_res = BlackBoxOptim.bboptimize(theta -> score_bthis(subsampled_id_list, data, theta, subsidy_type,
										              calibrated_delta = temp_calibrated_delta)[1];
		                                SearchRange = vcat(temp_search_domain[1:(param_dim-1)],(-100.0, 2000.0)),
										NumDimensions = param_dim,
										Method = :de_rand_1_bin,
										MaxSteps = num_steps_DE)
		#restore results
		#println("score: ", score_bthis(subsampled_id_list, data, m_res.archive_output.best_candidate, subsidy_type)[2])
		~, num_correct_ineq[iter,1], num_total_ineq[iter,1] = score_bthis(subsampled_id_list,
		                                                                    data,
																			m_res.archive_output.best_candidate,
																			subsidy_type)
		#println("the number of correct inequalities: ", num_correct_ineq[iter,1])
		#println("the total number of inequalities: ", num_total_ineq[iter,1])
		myests[iter,:] = m_res.archive_output.best_candidate
	end
	return num_correct_ineq, num_total_ineq, myests
end

#point estimate
function score_bthis_full_X(subsampled_id_list, data, theta, subsidy_type;
	                        calibrated_delta = 1)
	#target_theta = vcat(1,theta) # first parameter must be normalized to 1
    target_theta = vcat(theta,
	                    calibrated_delta) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data, target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end
function score_bthis_scale_X_only(subsampled_id_list, data, theta, subsidy_type;
	                        calibrated_delta = 1)
	#target_theta = vcat(1, theta[1:3], zeros(4),theta[4:5]) # first parameter must be normalized to 1
	target_theta = vcat(theta[1:4], zeros(4),
	                    theta[5], calibrated_delta) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data, target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end
function score_bthis_scale_and_scope_X_only(subsampled_id_list, data, theta, subsidy_type;
	                        calibrated_delta = 1)
	#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
	target_theta = vcat(zeros(4), theta,
	                    calibrated_delta) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data, target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end

function score_bthis_x45678_merger_cost(subsampled_id_list, data, theta, subsidy_type;
	                        calibrated_delta = 1)
	#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
	target_theta = vcat(zeros(3),
	                    theta[1:6],
	                    calibrated_delta) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data, target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end

# check behavior
temp_subsidy_type = "shared"
@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
                                                      score_bthis_full_X,
													  temp_subsidy_type,
													  temp_calibrated_delta = 1,
                                                      param_dim = 9,
													  num_its=1, size_of_subsample=60,
													  num_steps_DE=2)
@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
                                                      score_bthis_scale_X_only,
													  temp_subsidy_type,
													  temp_calibrated_delta = 1,
                                                      param_dim = 5,
													  num_its=1, size_of_subsample=60,
													  num_steps_DE=2)
@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
                                                      score_bthis_scale_and_scope_X_only,
													  temp_subsidy_type,
													  temp_calibrated_delta = 1,
                                                      param_dim = 5,
													  num_its=1, size_of_subsample=60,
													  num_steps_DE=2)
@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
                                                      score_bthis_x45678_merger_cost,
													  temp_subsidy_type,
													  temp_calibrated_delta = 1,
                                                      param_dim = 7,
													  num_its=1, size_of_subsample=106,
													  num_steps_DE=2)
#----------------------------#
temp_subsidy_type = "shared"
variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
function plot_score_single_variable(temp_subsidy_type;
	                                domain,
									data,
									variable_list,
									calibrated_delta_list,
	                                variable::String,
									file_name_variable::String,
									info_sum = temp_info_sum)
	#global data, variable_list
	single_variable_index = [1:1:length(variable_list);][variable_list .== variable][1]
	for kk in 1:length(calibrated_delta_list)
		calibrated_delta = calibrated_delta_list[kk]
		function score_bthis_only_target_x(subsampled_id_list, data,
			                               theta, subsidy_type)
			#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
			target_theta = vcat(zeros(single_variable_index-1),
			                    theta,
								zeros(9-single_variable_index),
								calibrated_delta) # first parameter must be normalized to 1
			score_res, total_num_ineq = score_b_est_data(subsampled_id_list,
			                            data,
			                            target_theta,
										subsidy_type,
										info_sum = temp_info_sum)
			res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
			println("score:$(res) \n" )
			return res, score_res, total_num_ineq
		end
		res_list = zeros(length(domain))
		for i in 1:length(domain)
			iter = domain[i]
			res_list[i], ~, total_num_ineq = score_bthis_only_target_x(data.firm_id,
			                             data, [iter],
										 temp_subsidy_type)
			global total_num_ineq
		end
		if kk == 1
			Plots.plot(title = "Objective value: subsidy type = $(temp_subsidy_type)",
			           legend = :bottomright,
					   ylim = [maximum(Int.(100000.0.-res_list))-50, maximum(Int.(100000.0.-res_list))+1],
			           xlabel = "$(variable)")
		end
		temp_accuracy = round(maximum(Int.(100000.0.-res_list))/total_num_ineq,digits=3)
		Plots.plot!(domain, Int.(100000.0.-res_list),
				   label = "score: its maximum of $(maximum(Int.(100000.0.-res_list))) out of $(total_num_ineq) inequalities ($(temp_accuracy)%),δ=$(calibrated_delta)")
		Plots.vline!(domain[-res_list.==maximum(-res_list)],
		             label = "maximum score: δ=$(calibrated_delta)",
					 color = "black",
					 linestyle = :dot)
	end
	savefig("julia_merger_figure/plot_score_only_$(file_name_variable)_$(temp_subsidy_type)_subsidy")
end

want_to_run = "not_run"
@time if want_to_run == "run"
	temp_subsidy_type = "shared"
	variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
	# economies of scale
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
	                           calibrated_delta_list = [0],
	                           variable = "β₁", file_name_variable = "x1")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
	                           calibrated_delta_list = [0],
	                           variable = "β₂", file_name_variable = "x2")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
	                           calibrated_delta_list = [0],
	                           variable = "β₃", file_name_variable = "x3")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = [0],
	                           variable = "β₄", file_name_variable = "x4")
	# economies of scope
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.5:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₅", file_name_variable = "x5")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.5:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₆", file_name_variable = "x6")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.5:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₇", file_name_variable = "x7")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₈", file_name_variable = "x8")
	# merger cost
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-1:10:400;],
							   calibrated_delta_list = [1,5,10],
	                           variable = "γ", file_name_variable = "gamma")
	temp_subsidy_type = "to_buyer"
	variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
	# economies of scale
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = [0],
	                           variable = "β₁", file_name_variable = "x1")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = [0],
	                           variable = "β₂", file_name_variable = "x2")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = [0],
	                           variable = "β₃", file_name_variable = "x3")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = [0],
	                           variable = "β₄", file_name_variable = "x4")
	# economies of scope
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₅", file_name_variable = "x5")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₆", file_name_variable = "x6")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₇", file_name_variable = "x7")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = [0],
	                           variable = "β₈", file_name_variable = "x8")
	# merger cost
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-1:5:400;],
							   calibrated_delta_list = [1,5,10],
	                           variable = "γ", file_name_variable = "gamma")
end
#730.659082 seconds (4.44 G allocations: 431.315 GiB, 7.14% gc time)


function plot_contour_score_two_variables(temp_subsidy_type;
	                                domain,
									data,
									variable_list,
									gamma_list,
									calibrated_delta_list,
	                                variable::String,
									file_name_variable::String,
									info_sum = temp_info_sum)
	#global data, variable_list, gamma_list
	single_variable_index = [1:1:length(variable_list);][variable_list .== variable][1]
	for kk in 1:length(calibrated_delta_list)
		calibrated_delta = calibrated_delta_list[kk]
		function score_bthis_only_target_x(subsampled_id_list,
			                               data,
			                               theta,
										   gg,
										   subsidy_type)
			#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
			target_theta = vcat(zeros(single_variable_index-1),
			                    theta,
								zeros(8-single_variable_index),
								gg,
								calibrated_delta) # first parameter must be normalized to 1
			score_res, total_num_ineq = score_b_est_data(subsampled_id_list,
										data,
										target_theta,
										subsidy_type,
										info_sum = temp_info_sum)
			res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
			println("score:$(res) \n" )
			return res, score_res, total_num_ineq
		end
		res_contor = zeros(length(domain),length(gamma_list))
		for i in 1:length(domain), gg in 1:length(gamma_list)
			iter = domain[i]
			gamma_iter = gamma_list[gg]
			res_contor[i,gg], ~, total_num_ineq = score_bthis_only_target_x(data.firm_id,
			                             data, [iter], [gamma_iter],
										 temp_subsidy_type)
			global total_num_ineq
		end
		if kk == 1
			Plots.contour(domain, gamma_list, Int.(100000.0.-res_contor)', fill = true,
			#Plots.plot(domain, gamma_list, Int.(100000.0.-res_contor)',st = [:contourf],
			              ylabel = "γ (merger cost)",
						  xlabel = "$variable",
						  title = "Objective value: (subsidy type)=($(temp_subsidy_type))")
		end
		#Plots.vline!([0], label="", linestyle = :dash)
	end
	savefig("julia_merger_figure/plot_contour_score_two_variables_merger_cost_$(file_name_variable)_$(temp_subsidy_type)_subsidy")
end


want_to_run = "not_run"
@time if want_to_run == "run"
	temp_calibrated_delta_list = [5]
	#gamma_list = [70:10.0:200;]
	#temp_domain = [100.0:50.0:1500;]
	temp_subsidy_type = "shared"
	variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
	gamma_list = [30:10.0:500;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-10.0:10.0:100;],
							   #domain = [-40.0:10:40;],
							   domain = [-20:20:200;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₂", file_name_variable = "x2")
	gamma_list = [60:10.0:500;]
	#gamma_list = [15:5.0:200;]
	#gamma_list = [19.52:0.02:19.62;]
	#gamma_list = [400:200.0:2000;]
	#temp_domain = [1100.0:50.0:1300;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-40:10:40;],
							   domain = [-10:2:10;],
							   #domain = [-9.91:0.01:-9.8;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₄", file_name_variable = "x4")
	#gamma_list = [20.0:2.0:200.0;]
	#gamma_list = [100:50.0:300;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-20:3:1;],
							   domain = [-40:10:0;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₅", file_name_variable = "x5")
	#gamma_list = [30.0:10.0:200.0;]
	#gamma_list = [50:10.0:130;]
	gamma_list = [80:10.0:500;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-10:2:0;],
							   domain = [-20:10:0;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₇", file_name_variable = "x7")
	#gamma_list = [40.0:10.0:200.0;]
	#gamma_list = [80:10.0:200;]
	#gamma_list = [50:50.0:300;]
	gamma_list = [70:2.0:120;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:5:15;],
							   domain = [7:0.3:19;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₈", file_name_variable = "x8")

	# change scale of merger cost
	#gamma_list = [0.4:0.2:1.5;]
	#gamma_list = [69:0.1:72.0;]
	#gamma_list = [60:5.0:120;]
	#gamma_list = [40:10.0:100;]
	gamma_list = [380:5.0:430;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:5:15;],
							   domain = [90:2:100;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₆", file_name_variable = "x6")
	gamma_list = [300:50.0:2000;]
	#gamma_list = [0.3:1.0:10.0;]
	#gamma_list = [30.0:2.0:200.0;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:10:40;],
							   domain = [10:20:200;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₁", file_name_variable = "x1")
	#gamma_list = [60.0:1:250;]
	gamma_list = [300:50.0:2000;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-20:10:40;],
							   domain = [-40:20:200;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₃", file_name_variable = "x3")
end
#3039.685852 seconds (11.68 G allocations: 1.109 TiB, 4.31% gc time)

#---------------------------#
# point estimate
#---------------------------#

size_of_fullsample = length(data.firm_id)-12
function point_estimate(subsidy_type;
	                    size_of_fullsample = 106,
	                    num_steps_DE_temp = 50,
	                    num_its_temp = 100,
						temp_temp_calibrated_delta = 1)
	# model 1
    @time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                  score_bthis_scale_X_only,
													  subsidy_type,
													  temp_calibrated_delta = temp_temp_calibrated_delta,
                                                      param_dim = 5,
													  num_its = num_its_temp,
													  size_of_subsample = size_of_fullsample,
													  num_steps_DE = num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_scale_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_scale_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_scale_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
	# model 2
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                      score_bthis_scale_and_scope_X_only,
														  subsidy_type,
														  temp_calibrated_delta = temp_temp_calibrated_delta,
	                                                      param_dim = 5,
														  num_its=num_its_temp,
														  size_of_subsample=size_of_fullsample,
														  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_scope_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_scope_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_scope_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
	# model 3
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                      score_bthis_full_X,
														  subsidy_type,
														  temp_calibrated_delta = temp_temp_calibrated_delta,
	                                                      param_dim = 9,
														  num_its=num_its_temp,
														  size_of_subsample=size_of_fullsample,
														  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_full_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_full_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_full_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
	# model 4
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                      score_bthis_x45678_merger_cost,
														  subsidy_type,
														  temp_calibrated_delta = temp_temp_calibrated_delta,
	                                                      param_dim = 6,
														  num_its=num_its_temp,
														  size_of_subsample=size_of_fullsample,
														  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_x45678_merger_cost.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_x45678_merger_cost.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_x45678_merger_cost.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
end

#----------------------------------------------------------#
# construct 95 percent Confidence interval via bootstrap
#----------------------------------------------------------#
function construct_CI(subsidy_type;
	                  num_its_bootstrap = 100,
	                  num_steps_DE_temp = 50,
	                  size_of_subsample_temp = 60,
					  temp_temp_calibrated_delta = 1)
    # model 1
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                  score_bthis_scale_X_only,
													  subsidy_type,
													  temp_calibrated_delta = temp_temp_calibrated_delta,
                                                      param_dim = 5,
													  num_its=num_its_bootstrap,
													  size_of_subsample=size_of_subsample_temp,
													  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_scale_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_scale_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_scale_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
	# model 2
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                      score_bthis_scale_and_scope_X_only,
														  subsidy_type,
														  temp_calibrated_delta = temp_temp_calibrated_delta,
	                                                      param_dim = 5,
														  num_its=num_its_bootstrap,
														  size_of_subsample=size_of_subsample_temp,
														  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_scope_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_scope_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_scope_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
	# model 3
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                      score_bthis_full_X,
														  subsidy_type,
														  temp_calibrated_delta = temp_temp_calibrated_delta,
	                                                      param_dim = 9,
														  num_its=num_its_bootstrap,
														  size_of_subsample=size_of_subsample_temp,
														  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_full_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_full_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_full_X_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
	# model 4
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                      score_bthis_x45678_merger_cost,
														  subsidy_type,
														  temp_calibrated_delta = temp_temp_calibrated_delta,
	                                                      param_dim = 6,
														  num_its=num_its_bootstrap,
														  size_of_subsample=size_of_subsample_temp,
														  num_steps_DE=num_steps_DE_temp)
	open("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_x45678_merger_cost.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_x45678_merger_cost.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_x45678_merger_cost.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
end

#=
export JULIA_NUM_THREADS=6
using Base.Threads
@show nthreads()
data = zeros(Int64,nthreads())
function test!(data)
    id = threadid()
    n = nthreads()
    nmax = 100
    nbun = div(nmax,n)
    ista = (id-1)*nbun + 1
    iend = ista + nbun -1
    println("$ista $iend at $id")
    a = 0
    for i=ista:iend
        a += i
    end

    data[id] = a
end
Threads.threading_run() do
    @show i = threadid()
    @show data[i] = test!(data)
end
=#
you_want_run = "not_run"
if you_want_run == "run"
	temp_subsidy_type = "shared"
	@time point_estimate(temp_subsidy_type,
	                     num_steps_DE_temp = 100,#50,
	                     num_its_temp = 50,#20,
						 #num_its_temp = 2,
						 size_of_fullsample = 106,
						 temp_temp_calibrated_delta = 5)
	#@time point_estimate(num_steps_DE_temp = 100, num_its_temp = 2)
	#1007.291976 seconds (5.74 G allocations: 557.036 GiB, 6.62% gc time)
	#@time point_estimate(num_steps_DE_temp = 100,  num_its_temp = 10)
	#3812.250026 seconds (20.37 G allocations: 2.009 TiB, 6.69% gc time)
	#@time point_estimate(num_steps_DE_temp = 50,  num_its_temp = 100)
	#151653.998867 seconds (160.68 G allocations: 15.258 TiB, 6.47% gc time)
	@time construct_CI(temp_subsidy_type,
	                   num_its_bootstrap = 200,
	                   #num_its_bootstrap = 2,
	                   num_steps_DE_temp = 100,
					   size_of_subsample_temp = 30,
					   temp_temp_calibrated_delta = 5)
	#34558.137605 seconds (162.34 G allocations: 15.437 TiB, 5.96% gc time)
	#@time point_estimate(num_steps_DE_temp = 50,  num_its_bootstrap = 200)
	#73597.442566 seconds (90.89 G allocations: 8.674 TiB, 11.04% gc time)
end

you_want_run = "not_run"
if you_want_run == "run"
	temp_subsidy_type = "to_buyer"
	@time point_estimate(temp_subsidy_type,
	                     num_steps_DE_temp = 100,
	                     #num_its_temp = 2,
						 num_its_temp = 50,
						 size_of_fullsample = 106)
	#@time point_estimate(num_steps_DE_temp = 100, num_its_temp = 1)
	#246.815583 seconds (1.28 G allocations: 124.312 GiB, 5.89% gc time)
	#@time point_estimate(num_steps_DE_temp = 100,  num_its_temp = 10)
	#3812.250026 seconds (20.37 G allocations: 2.009 TiB, 6.69% gc time)
	#DE100*iter50: 288975.892965 seconds (254.48 G allocations: 18.570 TiB, 3.70% gc time, 0.00% compilation time)
	@time construct_CI(temp_subsidy_type,
	                   num_its_bootstrap = 200,
	                   #num_its_bootstrap = 2,
	                   num_steps_DE_temp = 100,
					   size_of_subsample_temp = 30)
	#34558.137605 seconds (162.34 G allocations: 15.437 TiB, 5.96% gc time)
	#200boot*100DE: 328593.725486 seconds (290.70 G allocations: 21.262 TiB, 4.26% gc time, 0.00% compilation time)
end
#--------------------#
# read output (shared subsidy)
#--------------------#
# model 1
temp_subsidy_type = "shared"
size_of_subsample_temp = 30
function generate_score_table_model_1234(;temp_subsidy_type = "shared",
	                                     size_of_subsample_temp = 60)
	myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
	num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
	num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
	accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
	final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
	# model 2
	myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
	num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
	num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
	accuracy_scope_X_only = vec(num_correct_ineq_scope_X_only./num_total_ineq_scope_X_only)
	final_ests_point_scope_X_only = round.(myests_point_scope_X_only[findmax(accuracy_scope_X_only)[2],:],digits=2)
	# model 3
	myests_point_full_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
	num_correct_ineq_full_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
	num_total_ineq_full_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
	accuracy_full_X_only = vec(num_correct_ineq_full_X_only./num_total_ineq_full_X_only)
	final_ests_point_full_X_only = round.(myests_point_full_X_only[findmax(accuracy_full_X_only)[2],:],digits=2)
	# model 4
	myests_point_x45678_merger_cost = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)
	num_correct_ineq_x45678_merger_cost = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)
	num_total_ineq_x45678_merger_cost = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)
	accuracy_x45678_merger_cost = vec(num_correct_ineq_x45678_merger_cost./num_total_ineq_x45678_merger_cost)
	final_ests_point_x45678_merger_cost = round.(myests_point_x45678_merger_cost[findmax(accuracy_x45678_merger_cost)[2],:],digits=2)

	# CI
	#beta = theta[1:8]
	#gamma = theta[9] # coefficient on merger cost
	#delta = theta[10] # coefficient on subsidy indicator
	myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
	myests_CI_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
	myests_CI_full_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
	myests_CI_x45678_merger_cost = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)

	CI_scale_X_only_table = round.(hcat(
	    Statistics.quantile(myests_CI_scale_X_only[:,1], [0.025,0.975]),
		Statistics.quantile(myests_CI_scale_X_only[:,2], [0.025,0.975]),
		Statistics.quantile(myests_CI_scale_X_only[:,3], [0.025,0.975]),
		Statistics.quantile(myests_CI_scale_X_only[:,4], [0.025,0.975]),
		Statistics.quantile(-myests_CI_scale_X_only[:,5], [0.025,0.975])
		),digits=2)
	CI_scope_X_only_table = round.(hcat(
	    Statistics.quantile(myests_CI_scope_X_only[:,1], [0.025,0.975]),
		Statistics.quantile(myests_CI_scope_X_only[:,2], [0.025,0.975]),
		Statistics.quantile(myests_CI_scope_X_only[:,3], [0.025,0.975]),
		Statistics.quantile(myests_CI_scope_X_only[:,4], [0.025,0.975]),
		Statistics.quantile(-myests_CI_scope_X_only[:,5], [0.025,0.975])
		),digits=2)
	CI_full_X_only_table = round.(hcat(
	    Statistics.quantile(myests_CI_full_X_only[:,1], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,2], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,3], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,4], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,5], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,6], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,7], [0.025,0.975]),
		Statistics.quantile(myests_CI_full_X_only[:,8], [0.025,0.975]),
		Statistics.quantile(-myests_CI_full_X_only[:,9], [0.025,0.975])
		),digits=2)
	CI_x45678_merger_cost_table = round.(hcat(
	    Statistics.quantile(myests_CI_x45678_merger_cost[:,1], [0.025,0.975]),
		Statistics.quantile(myests_CI_x45678_merger_cost[:,2], [0.025,0.975]),
		Statistics.quantile(myests_CI_x45678_merger_cost[:,3], [0.025,0.975]),
		Statistics.quantile(myests_CI_x45678_merger_cost[:,4], [0.025,0.975]),
		Statistics.quantile(myests_CI_x45678_merger_cost[:,5], [0.025,0.975]),
		Statistics.quantile(-myests_CI_x45678_merger_cost[:,6], [0.025,0.975])
		),digits=2)
	total_ineq_all = Int64(num_total_ineq_x45678_merger_cost[findmax(accuracy_x45678_merger_cost)[2]])
	LaTeXTabulars.latex_tabular("julia_merger_table/score_results_temp_$(temp_subsidy_type)_subsidy.tex",
	              Tabular("@{\\extracolsep{5pt}}lccccc"),
	              [Rule(:top),
				   ["","","", "Value Function", "", ""],
	               ["","",
				   "Point Estimate", "Point Estimate", "Point Estimate", "Point Estimate"],
				   ["","",
				   "[95\\% CI Set Identified]", "[95\\% CI Set Identified]", "[95\\% CI Set Identified]", "[95\\% CI Set Identified]"],
				   Rule(:mid),

				   ["Measure of economies of scale", "", "", "", "", ""],
				   #beta_0
				   ["", "", "", "", ""],
				   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0", "+1", "+1", "+1", "+1"],
				   ["" , "" , "(S)", "(S)", "(S)", "(S)"],
				   #beta_1
				   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
				    final_ests_point_scale_X_only[1], "", final_ests_point_full_X_only[1], ""],
				   ["" , "" ,
				   "[$(CI_scale_X_only_table[1,1]), $(CI_scale_X_only_table[2,1])]",
				   "",
				   "[$(CI_full_X_only_table[1,1]), $(CI_full_X_only_table[2,1])]", ""],
				   #beta_2
	               ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
				    final_ests_point_scale_X_only[2], "", final_ests_point_full_X_only[2],
					""],
				   ["" , "" ,
				   "[$(CI_scale_X_only_table[1,2]), $(CI_scale_X_only_table[2,2])]",
				   "",
				   "[$(CI_full_X_only_table[1,2]), $(CI_full_X_only_table[2,2])]",
				   ""],
				   #beta_3
	               ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
				    final_ests_point_scale_X_only[3], "",
					 final_ests_point_full_X_only[3], ""],
				   ["" , "" ,
				   "[$(CI_scale_X_only_table[1,3]), $(CI_scale_X_only_table[2,3])]",
				   "",
				   "[$(CI_full_X_only_table[1,3]), $(CI_full_X_only_table[2,3])]", ""],
				   # beta_4
				   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
				    final_ests_point_scale_X_only[4], "", final_ests_point_full_X_only[4],
					 final_ests_point_x45678_merger_cost[1]],
				   ["" , "" ,
				   "[$(CI_scale_X_only_table[1,4]), $(CI_scale_X_only_table[2,4])]",
				   "",
				   "[$(CI_full_X_only_table[1,4]), $(CI_full_X_only_table[2,4])]",
				   "[$(CI_x45678_merger_cost_table[1,1]), $(CI_x45678_merger_cost_table[2,1])]"],
				   #beta_5
				   ["Measure of economies of scope", "", "", "", "", ""],
				   ["", "", "", "", "", ""],
	               ["share of liner\$_{b}\$ \$\\times\$ share of liner\$_{t}\$", L"\beta_5",
				    "", final_ests_point_scope_X_only[1], final_ests_point_full_X_only[5],
					 final_ests_point_x45678_merger_cost[2]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,1]), $(CI_scope_X_only_table[2,1])]",
				   "[$(CI_full_X_only_table[1,5]), $(CI_full_X_only_table[2,5])]",
				   "[$(CI_x45678_merger_cost_table[1,2]), $(CI_x45678_merger_cost_table[2,2])]"],
				   #beta_6
				   ["share of tramper\$_{b}\$ \$\\times\$ share of tramper\$_{t}\$", L"\beta_6",
				    "", final_ests_point_scope_X_only[2], final_ests_point_full_X_only[6],
					 final_ests_point_x45678_merger_cost[3]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,2]), $(CI_scope_X_only_table[2,2])]",
				   "[$(CI_full_X_only_table[1,6]), $(CI_full_X_only_table[2,6])]",
				   "[$(CI_x45678_merger_cost_table[1,3]), $(CI_x45678_merger_cost_table[2,3])]"],
				   #beta_7
	               ["share of special\$_{b}\$ \$\\times\$ share of special\$_{t}\$", L"\beta_7",
				    "", final_ests_point_scope_X_only[3], final_ests_point_full_X_only[7],
					 final_ests_point_x45678_merger_cost[4]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,3]), $(CI_scope_X_only_table[2,3])]",
				   "[$(CI_full_X_only_table[1,7]), $(CI_full_X_only_table[2,7])]",
				   "[$(CI_x45678_merger_cost_table[1,4]), $(CI_x45678_merger_cost_table[2,4])]"],
				   #beta_8
	               ["share of tanker\$_{b}\$ \$\\times\$ share of tanker\$_{t}\$", L"\beta_8",
				   "", final_ests_point_scope_X_only[4], final_ests_point_full_X_only[8],
				   final_ests_point_x45678_merger_cost[5]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,4]), $(CI_scope_X_only_table[2,4])]",
				   "[$(CI_full_X_only_table[1,8]), $(CI_full_X_only_table[2,8])]",
				   "[$(CI_x45678_merger_cost_table[1,5]), $(CI_x45678_merger_cost_table[2,5])]"],
				   #gamma merger cost (note the sign is reverse)
				   ["", "", "", "", "", ""],
				   ["", "", "", "", "", ""],
	               ["merger cost", "-\$\\gamma\$",
				   -final_ests_point_scale_X_only[5], -final_ests_point_scope_X_only[5],
				    -final_ests_point_full_X_only[9], -final_ests_point_x45678_merger_cost[6]],
				   ["" , "" ,
				   "[$(CI_scale_X_only_table[1,5]), $(CI_scale_X_only_table[2,5])]",
				   "[$(CI_scope_X_only_table[1,5]), $(CI_scope_X_only_table[2,5])]",
				   "[$(CI_full_X_only_table[1,9]), $(CI_full_X_only_table[2,9])]",
				   "[$(CI_x45678_merger_cost_table[1,6]), $(CI_x45678_merger_cost_table[2,6])]"],
				   #delta subsidy sensitivity
	               ["subsidy sensitivity (\$s^{\\text{shared}}\$)", L"\delta",
				   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
				   5, 5, 5, 5],
				   ["", "", "", "", "", ""],
	               Rule(),           # a nice \hline to make it ugly
				   ["\$\\sharp\$ Inequalities in Point Estimate" , "" ,
				    total_ineq_all,
					total_ineq_all,
					total_ineq_all,
					total_ineq_all
					],
				   ["\\% Inequalities" , "" ,
				   round(accuracy_scale_X_only[findmax(accuracy_scale_X_only)[2]],digits=3),
				   round(accuracy_scope_X_only[findmax(accuracy_scope_X_only)[2]],digits=3),
				   round(accuracy_full_X_only[findmax(accuracy_full_X_only)[2]],digits=3),
				   round(accuracy_x45678_merger_cost[findmax(accuracy_x45678_merger_cost)[2]],digits=3)],
	               Rule(:bottom)])
end
you_want_run = "not_run"
if you_want_run == "run"
    generate_score_table_model_1234(temp_subsidy_type = "shared",
	                                     size_of_subsample_temp = 30#60
										 )
    generate_score_table_model_1234(temp_subsidy_type = "to_buyer",
	                                    size_of_subsample_temp = 30)
end

#-----------------------------------------------------#
# find point-estimate LB of model 2
#-----------------------------------------------------#
function plot_point_estimate_histogram(temp_subsidy_type)
	myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)
	num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)
	num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)
	accuracy_scope_X_only = vec(num_correct_ineq_scope_X_only./num_total_ineq_scope_X_only)
	final_ests_point_scope_X_only = round.(myests_point_scope_X_only[findmax(accuracy_scope_X_only)[2],:],digits=2)
	best_accuracy_scope_X_only = round.(findmax(accuracy_scope_X_only)[1],digits=7)
	# second best
	temp_length = [1:1:length(accuracy_scope_X_only);]
	myests_point_scope_X_only_excluding_best = myests_point_scope_X_only[temp_length.!=findmax(accuracy_scope_X_only)[2],:]
	accuracy_scope_X_only_excluding_best = accuracy_scope_X_only[temp_length.!=findmax(accuracy_scope_X_only)[2]]
	second_best_ests_point_scope_X_only = round.(myests_point_scope_X_only_excluding_best[findmax(accuracy_scope_X_only_excluding_best)[2],:],digits=2)
	second_best_accuracy_scope_X_only = round.(findmax(accuracy_scope_X_only_excluding_best)[1],digits=7)
	# third best
	temp_length = [1:1:length(accuracy_scope_X_only_excluding_best);]
	myests_point_scope_X_only_excluding_best_and_secondbest = myests_point_scope_X_only_excluding_best[temp_length.!=findmax(accuracy_scope_X_only_excluding_best)[2],:]
	accuracy_scope_X_only_excluding_best_and_secondbest = accuracy_scope_X_only_excluding_best[temp_length.!=findmax(accuracy_scope_X_only)[2]]
	third_best_ests_point_scope_X_only = round.(myests_point_scope_X_only_excluding_best_and_secondbest[findmax(accuracy_scope_X_only_excluding_best_and_secondbest)[2],:],digits=2)
	third_best_accuracy_scope_X_only = round.(findmax(accuracy_scope_X_only_excluding_best_and_secondbest)[1],digits=7)

	Plots.histogram(myests_point_scope_X_only[:,1], bins=40,
	                ylabel = "count",
	                xlabel = "β₄",
					title = "Subsidy type ($temp_subsidy_type)",
	                label = "point-estimated β₄",
					xlim =[-200,200],
					alpha = 0.3)
	Plots.vline!([final_ests_point_scope_X_only[1]],
	             label = "β₄ (highest score:$(best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "black")
	Plots.vline!([second_best_ests_point_scope_X_only[1]],
	             label = "β₄ (second best score:$(second_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "red")
	Plots.vline!([third_best_ests_point_scope_X_only[1]],
	             label = "β₄ (third best score:$(third_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "blue")
	savefig("julia_merger_figure/histogram_estimated_parameters_beta_4_scope_X_only_$(temp_subsidy_type)_subsidy")

	Plots.histogram(myests_point_scope_X_only[:,2], bins=40,
	                ylabel = "count",
	                xlabel = "β₅",
					title = "Subsidy type ($temp_subsidy_type)",
	                label = "point-estimated β₅",
					xlim =[-200,200],
					alpha = 0.3)
	Plots.vline!([final_ests_point_scope_X_only[2]],
	             label = "β₅ (highest score:$(best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "black")
	Plots.vline!([second_best_ests_point_scope_X_only[2]],
	             label = "β₅ (second best score:$(second_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "red")
	Plots.vline!([third_best_ests_point_scope_X_only[2]],
	             label = "β₅ (third best score:$(third_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "blue")
	savefig("julia_merger_figure/histogram_estimated_parameters_beta_5_scope_X_only_$(temp_subsidy_type)_subsidy")
	Plots.histogram(myests_point_scope_X_only[:,3], bins=40,
	                ylabel = "count",
	                xlabel = "β₆",
					title = "Subsidy type ($temp_subsidy_type)",
	                label = "point-estimated β₆",
					xlim =[-200,200],
					alpha = 0.3)
	Plots.vline!([final_ests_point_scope_X_only[3]],
	             label = "β₆ (highest score:$(best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "black")
	Plots.vline!([second_best_ests_point_scope_X_only[3]],
	             label = "β₆ (second best score:$(second_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "red")
	Plots.vline!([third_best_ests_point_scope_X_only[3]],
	             label = "β₆ (third best score:$(third_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "blue")
	savefig("julia_merger_figure/histogram_estimated_parameters_beta_6_scope_X_only_$(temp_subsidy_type)_subsidy")
	Plots.histogram(myests_point_scope_X_only[:,4], bins=40,
	                ylabel = "count",
	                xlabel = "β₇",
					title = "Subsidy type ($temp_subsidy_type)",
	                label = "point-estimated β₇",
					xlim =[-200,200],
					alpha = 0.3)
	Plots.vline!([final_ests_point_scope_X_only[4]],
	             label = "β₇ (highest score:$(best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "black")
	Plots.vline!([second_best_ests_point_scope_X_only[4]],
	             label = "β₇ (second best score:$(second_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "red")
	Plots.vline!([third_best_ests_point_scope_X_only[4]],
	             label = "β₇ (third best score:$(third_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "blue")
	savefig("julia_merger_figure/histogram_estimated_parameters_beta_7_scope_X_only_$(temp_subsidy_type)_subsidy")
	Plots.histogram(myests_point_scope_X_only[:,5], bins=40,
	                ylabel = "count",
	                xlabel = "β₈",
					title = "Subsidy type ($temp_subsidy_type)",
	                label = "point-estimated β₈",
					xlim =[-200,200],
					alpha = 0.3)
	Plots.vline!([final_ests_point_scope_X_only[5]],
	             label = "β₈ (highest score:$(best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "black")
	Plots.vline!([second_best_ests_point_scope_X_only[5]],
	             label = "β₈ (second best score:$(second_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "red")
	Plots.vline!([third_best_ests_point_scope_X_only[5]],
	             label = "β₈ (third best score:$(third_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "blue")
	savefig("julia_merger_figure/histogram_estimated_parameters_beta_8_scope_X_only_$(temp_subsidy_type)_subsidy")
	Plots.histogram(-myests_point_scope_X_only[:,6], bins=40,
	                ylabel = "count",
	                xlabel = "γ (merger cost)",
					title = "Subsidy type ($temp_subsidy_type)",
	                label = "point-estimated γ",
					xlim =[-200,200],
					alpha = 0.3)
	Plots.vline!([-final_ests_point_scope_X_only[6]],
	             label = "γ (highest score:$(best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "black")
	Plots.vline!([-second_best_ests_point_scope_X_only[6]],
	             label = "γ (second best score:$(second_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "red")
	Plots.vline!([-third_best_ests_point_scope_X_only[6]],
	             label = "γ (third best score:$(third_best_accuracy_scope_X_only))",
				 linestyle = :dash,
	             color = "blue")
	savefig("julia_merger_figure/histogram_estimated_parameters_gamma_scope_X_only_$(temp_subsidy_type)_subsidy")
end
you_want_run = "not_run"
if you_want_run == "run"
    @time plot_point_estimate_histogram("shared")
	#0.201591 seconds (261.57 k allocations: 10.609 MiB)
    #plot_point_estimate_histogram("to_buyer")
end

#----------------------------------#
# two_variables
#----------------------------------#
# two_variables
temp_subsidy_type = "shared"
variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
function point_estimate_two_variables(subsidy_type;
	                    data,
						variable_list,
	                    size_of_fullsample = 106,
	                    num_steps_DE_temp = 50,
	                    num_its_temp = 100,
						calibrated_delta_list = [1],
						variable::String,
						file_name_variable::String,
						info_sum = temp_info_sum)
	#global data, variable_list, gamma_list
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
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
			score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data,
										target_theta, subsidy_type)
			res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
			#println("score:$(res) \n" )
			return res, score_res, total_num_ineq
		end
		@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
		                                                      score_bthis_only_target_x,
															  temp_subsidy_type,
															  temp_calibrated_delta = calibrated_delta_kk,
		                                                      param_dim = 2,
															  num_its=num_its_temp,
															  size_of_subsample=size_of_fullsample,
															  num_steps_DE=num_steps_DE_temp)
		open("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt", "w") do io
			DelimitedFiles.writedlm(io, myests,",")
		end
		open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt", "w") do io
			DelimitedFiles.writedlm(io, num_total_ineq,",")
		end
		open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt", "w") do io
			DelimitedFiles.writedlm(io, num_correct_ineq,",")
		end
	end
end

function construct_CI_two_variables(subsidy_type;
					  data,
					  variable_list,
	                  num_its_bootstrap = 100,
	                  num_steps_DE_temp = 50,
	                  size_of_subsample_temp = 60,
					  calibrated_delta_list = [1],
					  variable::String,
					  file_name_variable::String,
					  info_sum = temp_info_sum)
	#global data
	#global variable_list, gamma_list
	liner_sum = info_sum["liner_sum"]
	special_sum = info_sum["special_sum"]
	tanker_sum = info_sum["tanker_sum"]
	tramper_sum = info_sum["tramper_sum"]
	single_variable_index = [1:1:length(variable_list);][variable_list .== variable][1]
	for kk = 1:length(calibrated_delta_list)
		calibrated_delta_kk = calibrated_delta_list[kk]
		function score_bthis_only_target_x(subsampled_id_list, data,
										   theta, subsidy_type;
										   calibrated_delta = calibrated_delta_kk)
			#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
			target_theta = vcat(zeros(single_variable_index-1),
								theta[1],
								zeros(8-single_variable_index),
								theta[2],
								calibrated_delta) # first parameter must be normalized to 1
			score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data,
										target_theta, subsidy_type)
			res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
			#println("score:$(res) \n" )
			return res, score_res, total_num_ineq
		end
		@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
		                                                      score_bthis_only_target_x,
															  temp_subsidy_type,
															  temp_calibrated_delta = calibrated_delta_kk,
		                                                      param_dim = 2,
															  num_its=num_its_bootstrap,
															  size_of_subsample=size_of_subsample_temp,
															  num_steps_DE=num_steps_DE_temp)
		open("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt", "w") do io
			DelimitedFiles.writedlm(io, myests,",")
		end
		open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt", "w") do io
			DelimitedFiles.writedlm(io, num_total_ineq,",")
		end
		open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt", "w") do io
			DelimitedFiles.writedlm(io, num_correct_ineq,",")
		end
	end
end
#-----------------------#
# estimate shared case
#-----------------------#
temp_subsidy_type = "shared"
variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
file_name_variable_list = ["x1","x2","x3","x4","x5","x6","x7","x8"]
you_want_run = "not_run"
Threads.nthreads()
JULIA_NUM_THREADS=8
temp_calibrated_delta_list = [5]
if you_want_run == "run"
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
	    @time point_estimate_two_variables(temp_subsidy_type;
							data = data,
	    					variable_list = variable_list,
		                    size_of_fullsample = 106,
		                    num_steps_DE_temp = 100,
		                    num_its_temp = 50,
							calibrated_delta_list = temp_calibrated_delta_list,
							variable = temp_target_variable,
							file_name_variable = temp_file_name,
							info_sum = temp_info_sum)
	end
	#300*8model
	#3343.760520 seconds (34.51 G allocations: 3.277 TiB, 67.25% gc time)
	#5000*8model
	#87096.960875 seconds (320.63 G allocations: 30.447 TiB, 33.96% gc time)
	#2500*8models
	#44417.409025 seconds (285.45 G allocations: 20.823 TiB, 32.66% gc time, 0.00% compilation time)
	#1000*8models
	#6474.324370 seconds (112.92 G allocations: 8.237 TiB, 30.34% gc time)
	#100DE*50step*8models
	#76386.576026 seconds (502.18 G allocations: 36.632 TiB, 34.90% gc time, 0.00% compilation time)
end

if you_want_run == "run"
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
		construct_CI_two_variables(temp_subsidy_type;
							  data = data,
    	    				  variable_list = variable_list,
			                  num_its_bootstrap = 200,
			                  num_steps_DE_temp = 100,
			                  size_of_subsample_temp = 30,
							  calibrated_delta_list = temp_calibrated_delta_list,
							  variable = temp_target_variable,
							  file_name_variable = temp_file_name,
  							  info_sum = temp_info_sum)
	end
	#500*8model
	#668.661895 seconds (13.49 G allocations: 1.288 TiB, 35.19% gc time)
	#10000*8
	#13125.682549 seconds (180.28 G allocations: 17.204 TiB, 23.98% gc time)
	#140238.952268 seconds (576.93 G allocations: 41.820 TiB, 24.50% gc time, 0.00% compilation time)
	# 200boot*50DE on 210430
	#52134.231030 seconds (323.61 G allocations: 23.661 TiB, 20.94% gc time, 0.00% compilation time)
	#20015.040089 seconds (325.65 G allocations: 23.805 TiB, 27.57% gc time, 0.00% compilation time)
	#100DE*50step*8models
	#107468.343673 seconds (573.67 G allocations: 41.936 TiB, 17.23% gc time, 0.00% compilation time)
end
temp_subsidy_type = "to_buyer"
variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
file_name_variable_list = ["x1","x2","x3","x4","x5","x6","x7","x8"]
you_want_run = "not_run"
Threads.nthreads()
JULIA_NUM_THREADS=8
temp_calibrated_delta_list = [5]
if you_want_run == "run"
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
	    @time point_estimate_two_variables(temp_subsidy_type;
							data = data,
	    					variable_list = variable_list,
		                    size_of_fullsample = 106,
		                    num_steps_DE_temp = 100,
		                    num_its_temp = 50,
							calibrated_delta_list = temp_calibrated_delta_list,
							variable = temp_target_variable,
							file_name_variable = temp_file_name,
							info_sum = temp_info_sum)
	end
	#300*8model
	#3343.760520 seconds (34.51 G allocations: 3.277 TiB, 67.25% gc time)
	#5000*8model
	#87096.960875 seconds (320.63 G allocations: 30.447 TiB, 33.96% gc time)
	#2500*8models
	#44417.409025 seconds (285.45 G allocations: 20.823 TiB, 32.66% gc time, 0.00% compilation time)
	#1000*8models
	#6474.324370 seconds (112.92 G allocations: 8.237 TiB, 30.34% gc time)
	#100DE*50step*8models
	#76386.576026 seconds (502.18 G allocations: 36.632 TiB, 34.90% gc time, 0.00% compilation time)
end

if you_want_run == "run"
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
		construct_CI_two_variables(temp_subsidy_type;
							  data = data,
    	    				  variable_list = variable_list,
			                  num_its_bootstrap = 200,
			                  num_steps_DE_temp = 100,
			                  size_of_subsample_temp = 30,
							  calibrated_delta_list = temp_calibrated_delta_list,
							  variable = temp_target_variable,
							  file_name_variable = temp_file_name,
  							  info_sum = temp_info_sum)
	end
	#500*8model
	#668.661895 seconds (13.49 G allocations: 1.288 TiB, 35.19% gc time)
	#10000*8
	#13125.682549 seconds (180.28 G allocations: 17.204 TiB, 23.98% gc time)
	#140238.952268 seconds (576.93 G allocations: 41.820 TiB, 24.50% gc time, 0.00% compilation time)
	# 200boot*50DE on 210430
	#52134.231030 seconds (323.61 G allocations: 23.661 TiB, 20.94% gc time, 0.00% compilation time)
	#20015.040089 seconds (325.65 G allocations: 23.805 TiB, 27.57% gc time, 0.00% compilation time)
	#100DE*50step*8models
	#107468.343673 seconds (573.67 G allocations: 41.936 TiB, 17.23% gc time, 0.00% compilation time)
end

# check behavior
temp_subsidy_type = "shared"
file_name_variable = "x1"
size_of_subsample_temp = 30
myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)

model_all_length = length(file_name_variable_list)
myests_point_all = zeros(model_all_length,
                         size(myests_point_scale_X_only)[1],
						 size(myests_point_scale_X_only)[2])
num_correct_ineq_all = zeros(model_all_length,
                         size(num_correct_ineq_scale_X_only)[1],
						 size(num_correct_ineq_scale_X_only)[2])
num_total_ineq_all = zeros(model_all_length,
                         size(num_total_ineq_scale_X_only)[1],
						 size(num_total_ineq_scale_X_only)[2])
accuracy_all = zeros(model_all_length,
                         size(accuracy_scale_X_only)[1])
final_ests_point_all = zeros(model_all_length,
                         size(final_ests_point_scale_X_only)[1])
myests_CI_all = zeros(model_all_length,
                         size(myests_CI_scale_X_only)[1],
						 size(myests_CI_scale_X_only)[2])
for ii = 1:length(file_name_variable_list)
	temp_subsidy_type = "shared"
	@show temp_file_name = file_name_variable_list[ii]
	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(accuracy_all[ii,:])[2], :],
	                                    digits=2)
	myests_CI_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
end
maximum(num_correct_ineq_all[:,:,1],dims = 2)
final_ests_point_all
myests_CI_all[1,:,2]
# CI
#beta = theta[1:8]
#gamma = theta[9] # coefficient on merger cost
#delta = theta[10] # coefficient on subsidy indicator
CI_all_table = zeros(model_all_length,2,2)
for ii = 1:length(file_name_variable_list)
	CI_all_table[ii,:,:] = round.(hcat(
		Statistics.quantile(myests_CI_all[ii,:,1], [0.025,0.975]),
		Statistics.quantile(-myests_CI_all[ii,:,2], [0.025,0.975])
		),digits=1)
end
CI_all_table[1,:,:]
maximum(num_correct_ineq_all[1,:,:])
total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy.tex",
			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			  [Rule(:top),
			   ["","","", "", "", "", "", "", "", ""],
			   ["","",
			   "Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   ["","",
			   "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]",
			    "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]"],
			   Rule(:mid),
			   ["Economies of scale", "", "", "", "", "", "", ""],
			   #beta_0
			   ["", "", "", "", "", "", "", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
			    "(S)", "(S)", "(S)", "(S)"],
			   #beta_1
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
				final_ests_point_all[1,1], "", "", "",
				 "", "", "", ""],
			   ["" , "" ,
			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
			   "", "", "",
			   "", "", "", ""],
			   #beta_2
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
				"", final_ests_point_all[2,1], "", "",
				"", "", "", ""],
			   ["" , "" ,
			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
			   "", "", "", ""],
			   #beta_3
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			   "", "", final_ests_point_all[3,1], "",
			   "", "", "", ""],
			   ["" , "" ,
			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
			   "", "", "", ""],
			   # beta_4
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			   "", "", "", final_ests_point_all[4,1],
			   "", "", "", ""],
			   ["" , "" ,
			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
			   "", "", "", ""],
			   #beta_5
			   ["Economies of scope", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
			   "", "", "", "",
			   final_ests_point_all[5,1], "", "", ""],
			   ["" , "" ,
			   "", "", "", "",
			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
			   #beta_6
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
			   "", "", "", "",
			   "", final_ests_point_all[6,1], "", ""],
			   ["" , "" ,
			   "", "", "", "",
			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
			   #beta_7
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
			   "", "", "", "",
			   "", "", final_ests_point_all[7,1], ""],
			   ["" , "" ,
			   "", "", "", "",
			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
			   #beta_8
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
			   "", "", "", "",
			   "", "", "", final_ests_point_all[8,1]],
			   ["" , "" ,
			   "", "", "", "",
			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
			   #gamma merger cost (note the sign is reverse)
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["merger cost", "-\$\\gamma\$",
			   -final_ests_point_all[1,2], -final_ests_point_all[2,2],
			   -final_ests_point_all[3,2], -final_ests_point_all[4,2],
			   -final_ests_point_all[5,2], -final_ests_point_all[6,2],
			   -final_ests_point_all[7,2], -final_ests_point_all[8,2]],
			   ["" , "" ,
			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
			   #delta subsidy sensitivity
			   ["subsidy sensitivity", L"\delta",
			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
			   ["", "", "", "", "", "", "", "", "", ""],
			   Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all],
				#=["\\Correct Ineq" , "" ,
 			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
 			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
 			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
 			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
 			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
 			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
 			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
 			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
			   =#
			   ["\\% Inequalities" , "" ,
			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   Rule(:bottom)])

relative_coefficients = round.([final_ests_point_all[1,1]/abs(final_ests_point_all[1,2]),
			   final_ests_point_all[2,1]/abs(final_ests_point_all[2,2]),
			   final_ests_point_all[3,1]/abs(final_ests_point_all[3,2]),
			   final_ests_point_all[4,1]/abs(final_ests_point_all[4,2]),
			   final_ests_point_all[5,1]/abs(final_ests_point_all[5,2]),
			   final_ests_point_all[6,1]/abs(final_ests_point_all[6,2]),
			   final_ests_point_all[7,1]/abs(final_ests_point_all[7,2]),
			   final_ests_point_all[8,1]/abs(final_ests_point_all[8,2])],digits = 3)
LaTeXTabulars.latex_tabular("julia_merger_table/ratio_score_results_two_variables_$(temp_subsidy_type)_subsidy.tex",
			  Tabular("@{\\extracolsep{5pt}}lcccccccc"),
			  [Rule(:top),
			   ["","\$\\beta_1/|\\gamma|\$","\$\\beta_2/|\\gamma|\$",
			    "\$\\beta_3/|\\gamma|\$", "\$\\beta_4/|\\gamma|\$",
				"\$\\beta_5/|\\gamma|\$", "\$\\beta_6/|\\gamma|\$",
				"\$\\beta_7/|\\gamma|\$", "\$\\beta_8/|\\gamma|\$"],
			   Rule(:mid),
			   ["", "", "", "", "", "", "", "",""],
			   vcat("",
			   relative_coefficients),
			   ["", "", "", "", "", "", "", "",""],
			   #delta subsidy sensitivity
			   vcat("relative level",
			   round.(relative_coefficients./minimum(abs.(relative_coefficients)),digits = 2)),
			   ["", "", "", "", "", "", "", "",""],
			   Rule(:bottom)])
for ii = 1:length(file_name_variable_list)
	temp_subsidy_type = "to_buyer"
	@show temp_file_name = file_name_variable_list[ii]
	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(accuracy_all[ii,:])[2], :],
	                                    digits=2)
	myests_CI_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
end
maximum(num_correct_ineq_all[:,:,1],dims = 2)
final_ests_point_all
myests_CI_all[1,:,2]
# CI
#beta = theta[1:8]
#gamma = theta[9] # coefficient on merger cost
#delta = theta[10] # coefficient on subsidy indicator
CI_all_table = zeros(model_all_length,2,2)
for ii = 1:length(file_name_variable_list)
	CI_all_table[ii,:,:] = round.(hcat(
		Statistics.quantile(myests_CI_all[ii,:,1], [0.025,0.975]),
		Statistics.quantile(-myests_CI_all[ii,:,2], [0.025,0.975])
		),digits=1)
end
CI_all_table[1,:,:]
maximum(num_correct_ineq_all[1,:,:])
total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy.tex",
			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			  [Rule(:top),
			   ["","","", "", "", "", "", "", "", ""],
			   ["","",
			   "Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   ["","",
			   "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]",
			    "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]"],
			   Rule(:mid),
			   ["Economies of scale", "", "", "", "", "", "", ""],
			   #beta_0
			   ["", "", "", "", "", "", "", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
			    "(S)", "(S)", "(S)", "(S)"],
			   #beta_1
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
				final_ests_point_all[1,1], "", "", "",
				 "", "", "", ""],
			   ["" , "" ,
			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
			   "", "", "",
			   "", "", "", ""],
			   #beta_2
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
				"", final_ests_point_all[2,1], "", "",
				"", "", "", ""],
			   ["" , "" ,
			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
			   "", "", "", ""],
			   #beta_3
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			   "", "", final_ests_point_all[3,1], "",
			   "", "", "", ""],
			   ["" , "" ,
			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
			   "", "", "", ""],
			   # beta_4
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			   "", "", "", final_ests_point_all[4,1],
			   "", "", "", ""],
			   ["" , "" ,
			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
			   "", "", "", ""],
			   #beta_5
			   ["Economies of scope", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
			   "", "", "", "",
			   final_ests_point_all[5,1], "", "", ""],
			   ["" , "" ,
			   "", "", "", "",
			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
			   #beta_6
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
			   "", "", "", "",
			   "", final_ests_point_all[6,1], "", ""],
			   ["" , "" ,
			   "", "", "", "",
			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
			   #beta_7
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
			   "", "", "", "",
			   "", "", final_ests_point_all[7,1], ""],
			   ["" , "" ,
			   "", "", "", "",
			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
			   #beta_8
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
			   "", "", "", "",
			   "", "", "", final_ests_point_all[8,1]],
			   ["" , "" ,
			   "", "", "", "",
			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
			   #gamma merger cost (note the sign is reverse)
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["merger cost", "-\$\\gamma\$",
			   -final_ests_point_all[1,2], -final_ests_point_all[2,2],
			   -final_ests_point_all[3,2], -final_ests_point_all[4,2],
			   -final_ests_point_all[5,2], -final_ests_point_all[6,2],
			   -final_ests_point_all[7,2], -final_ests_point_all[8,2]],
			   ["" , "" ,
			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
			   #delta subsidy sensitivity
			   ["subsidy sensitivity", L"\delta",
			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
			   ["", "", "", "", "", "", "", "", "", ""],
			   Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all],
				#=["\\Correct Ineq" , "" ,
 			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
 			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
 			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
 			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
 			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
 			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
 			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
 			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
			   =#
			   ["\\% Inequalities" , "" ,
			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   Rule(:bottom)])

relative_coefficients = round.([final_ests_point_all[1,1]/abs(final_ests_point_all[1,2]),
			   final_ests_point_all[2,1]/abs(final_ests_point_all[2,2]),
			   final_ests_point_all[3,1]/abs(final_ests_point_all[3,2]),
			   final_ests_point_all[4,1]/abs(final_ests_point_all[4,2]),
			   final_ests_point_all[5,1]/abs(final_ests_point_all[5,2]),
			   final_ests_point_all[6,1]/abs(final_ests_point_all[6,2]),
			   final_ests_point_all[7,1]/abs(final_ests_point_all[7,2]),
			   final_ests_point_all[8,1]/abs(final_ests_point_all[8,2])],digits = 3)
LaTeXTabulars.latex_tabular("julia_merger_table/ratio_score_results_two_variables_$(temp_subsidy_type)_subsidy.tex",
			  Tabular("@{\\extracolsep{5pt}}lcccccccc"),
			  [Rule(:top),
			   ["","\$\\beta_1/|\\gamma|\$","\$\\beta_2/|\\gamma|\$",
			    "\$\\beta_3/|\\gamma|\$", "\$\\beta_4/|\\gamma|\$",
				"\$\\beta_5/|\\gamma|\$", "\$\\beta_6/|\\gamma|\$",
				"\$\\beta_7/|\\gamma|\$", "\$\\beta_8/|\\gamma|\$"],
			   Rule(:mid),
			   ["", "", "", "", "", "", "", "",""],
			   vcat("",
			   relative_coefficients),
			   ["", "", "", "", "", "", "", "",""],
			   #delta subsidy sensitivity
			   vcat("relative level",
			   round.(relative_coefficients./minimum(abs.(relative_coefficients)),digits = 2)),
			   ["", "", "", "", "", "", "", "",""],
			   Rule(:bottom)])

#-----------------------#
# estimate only main firm
#-----------------------#
temp_subsidy_type = "shared"
variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
file_name_variable_list = ["x1","x2","x3","x4","x5","x6","x7","x8"]
you_want_run = "not_run"
Threads.nthreads()
JULIA_NUM_THREADS=8
temp_calibrated_delta_list = [1000]
if you_want_run == "run"
	# choose only main firms
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
	    @time point_estimate_two_variables(temp_subsidy_type;
		                    data = data,
		                    variable_list = variable_list,
		                    size_of_fullsample = 0,
							num_steps_DE_temp = 200,#50,
		                    num_its_temp = 200,#1000,
							calibrated_delta_list = temp_calibrated_delta_list,
							variable = temp_target_variable,
							file_name_variable = temp_file_name,
							info_sum = temp_info_sum)
	end
	#300*8model
	#228.085250 seconds (4.92 G allocations: 364.622 GiB, 32.77% gc time)
	#50*1000*8model
	#35169.052760 seconds (49.28 G allocations: 3.564 TiB, 10.58% gc time)
	#200step*1000iter*8model
	#7607.351116 seconds (148.67 G allocations: 10.799 TiB, 30.05% gc time, 0.00% compilation time)
end
temp_subsidy_type = "to_buyer"
if you_want_run == "run"
	# choose only main firms
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
	    @time point_estimate_two_variables(temp_subsidy_type;
		                    data = data,
		                    variable_list = variable_list,
		                    size_of_fullsample = 0,
							num_steps_DE_temp = 200,#50,
		                    num_its_temp = 200,#1000,
							calibrated_delta_list = temp_calibrated_delta_list,
							variable = temp_target_variable,
							file_name_variable = temp_file_name,
							info_sum = temp_info_sum)
	end
	#300*8model
	#228.085250 seconds (4.92 G allocations: 364.622 GiB, 32.77% gc time)
	#50*1000*8model
	#35169.052760 seconds (49.28 G allocations: 3.564 TiB, 10.58% gc time)
	#200step*1000iter*8model
	#7607.351116 seconds (148.67 G allocations: 10.799 TiB, 30.05% gc time, 0.00% compilation time)
end

# check behavior
temp_subsidy_type = "shared"
file_name_variable = "x1"
size_of_subsample_temp = 0
myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)

model_all_length = length(file_name_variable_list)
myests_point_all = zeros(model_all_length,
                         size(myests_point_scale_X_only)[1],
						 size(myests_point_scale_X_only)[2])
num_correct_ineq_all = zeros(model_all_length,
                         size(num_correct_ineq_scale_X_only)[1],
						 size(num_correct_ineq_scale_X_only)[2])
num_total_ineq_all = zeros(model_all_length,
                         size(num_total_ineq_scale_X_only)[1],
						 size(num_total_ineq_scale_X_only)[2])
accuracy_all = zeros(model_all_length,
                         size(accuracy_scale_X_only)[1])
final_ests_point_all = zeros(model_all_length,
                         size(final_ests_point_scale_X_only)[1])
for ii = 1:length(file_name_variable_list)
	temp_subsidy_type = "to_buyer"
	@show temp_file_name = file_name_variable_list[ii]
	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(num_correct_ineq_all[ii,:,1])[2], :],
	                                    digits=1)
end

CI_all_table = zeros(model_all_length,2,2)
for ii = 1:length(file_name_variable_list)
	temp_myests_point_all = myests_point_all[ii,num_correct_ineq_all[ii,:,1].==Int(findmax(num_correct_ineq_all[ii,:,1])[1]),:]
	CI_all_table[ii,:,:] = round.(hcat(
		Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
		Statistics.quantile(-temp_myests_point_all[:,2], [0.0,1.0])
		),digits=1)
end

total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy_main_firms_only.tex",
			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			  [Rule(:top),
			   ["","","", "", "", "", "", "", "", ""],
			   #["","",
			   #"Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   ["","",
			   "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]",
			    "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]"],
			   Rule(:mid),
			   ["Economies of scale", "", "", "", "", "", "", ""],
			   #beta_0
			   ["", "", "", "", "", "", "", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
			    "(S)", "(S)", "(S)", "(S)"],
			   #beta_1
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
			   #final_ests_point_all[1,1], "", "", "",
			   # "", "", "", ""],
			   #["" , "" ,
			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
			   "", "", "",
			   "", "", "", ""],
			   #beta_2
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
				#"", final_ests_point_all[2,1], "", "",
				#"", "", "", ""],
			   #["" , "" ,
			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
			   "", "", "", ""],
			   #beta_3
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			   #"", "", final_ests_point_all[3,1], "",
			   #"", "", "", ""],
			   #["" , "" ,
			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
			   "", "", "", ""],
			   # beta_4
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			   #"", "", "", final_ests_point_all[4,1],
			   #"", "", "", ""],
			   #["" , "" ,
			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
			   "", "", "", ""],
			   #beta_5
			   ["Economies of scope", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
			   #"", "", "", "",
			   #final_ests_point_all[5,1], "", "", ""],
			   #["" , "" ,
			   "", "", "", "",
			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
			   #beta_6
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
			   #"", "", "", "",
			   #"", final_ests_point_all[6,1], "", ""],
			   #["" , "" ,
			   "", "", "", "",
			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
			   #beta_7
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
			   #"", "", "", "",
			   #"", "", final_ests_point_all[7,1], ""],
			   #["" , "" ,
			   "", "", "", "",
			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
			   #beta_8
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
			   #"", "", "", "",
			   #"", "", "", final_ests_point_all[8,1]],
			   #["" , "" ,
			   "", "", "", "",
			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
			   #gamma merger cost (note the sign is reverse)
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["merger cost", "-\$\\gamma\$",
			   #-final_ests_point_all[1,2], -final_ests_point_all[2,2],
			   #-final_ests_point_all[3,2], -final_ests_point_all[4,2],
			   #-final_ests_point_all[5,2], -final_ests_point_all[6,2],
			   #-final_ests_point_all[7,2], -final_ests_point_all[8,2]],
			   #["" , "" ,
			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
			   #delta subsidy sensitivity
			   ["subsidy sensitivity", L"\delta",
			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
			   ["", "", "", "", "", "", "", "", "", ""],
			   Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all],
				#=["\\Correct Ineq" , "" ,
 			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
 			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
 			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
 			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
 			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
 			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
 			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
 			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
			   =#
			   ["\\% Inequalities" , "" ,
			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   Rule(:bottom)])


			   # check behavior
			   temp_subsidy_type = "shared"
			   file_name_variable = "x1"
			   size_of_subsample_temp = 0
			   myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
			   num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
			   num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
			   accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
			   final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
			   myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)

			   model_all_length = length(file_name_variable_list)
			   myests_point_all = zeros(model_all_length,
			                            size(myests_point_scale_X_only)[1],
			   						 size(myests_point_scale_X_only)[2])
			   num_correct_ineq_all = zeros(model_all_length,
			                            size(num_correct_ineq_scale_X_only)[1],
			   						 size(num_correct_ineq_scale_X_only)[2])
			   num_total_ineq_all = zeros(model_all_length,
			                            size(num_total_ineq_scale_X_only)[1],
			   						 size(num_total_ineq_scale_X_only)[2])
			   accuracy_all = zeros(model_all_length,
			                            size(accuracy_scale_X_only)[1])
			   final_ests_point_all = zeros(model_all_length,
			                            size(final_ests_point_scale_X_only)[1])
			   for ii = 1:length(file_name_variable_list)
			   	temp_subsidy_type = "shared"
			   	@show temp_file_name = file_name_variable_list[ii]
			   	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
			   	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
			   	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
			   	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
			   	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(num_correct_ineq_all[ii,:,1])[2], :],
			   	                                    digits=1)
			   end

			   CI_all_table = zeros(model_all_length,2,2)
			   for ii = 1:length(file_name_variable_list)
			   	temp_myests_point_all = myests_point_all[ii,num_correct_ineq_all[ii,:,1].==Int(findmax(num_correct_ineq_all[ii,:,1])[1]),:]
			   	CI_all_table[ii,:,:] = round.(hcat(
			   		Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
			   		Statistics.quantile(-temp_myests_point_all[:,2], [0.0,1.0])
			   		),digits=1)
			   end

			   total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
			   LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy_main_firms_only.tex",
			   			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			   			  [Rule(:top),
			   			   ["","","", "", "", "", "", "", "", ""],
			   			   #["","",
			   			   #"Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   			   ["","",
			   			   "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]",
			   			    "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]"],
			   			   Rule(:mid),
			   			   ["Economies of scale", "", "", "", "", "", "", ""],
			   			   #beta_0
			   			   ["", "", "", "", "", "", "", "", ""],
			   			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
			   			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
			   			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
			   			    "(S)", "(S)", "(S)", "(S)"],
			   			   #beta_1
			   			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
			   			   #final_ests_point_all[1,1], "", "", "",
			   			   # "", "", "", ""],
			   			   #["" , "" ,
			   			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
			   			   "", "", "",
			   			   "", "", "", ""],
			   			   #beta_2
			   			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
			   				#"", final_ests_point_all[2,1], "", "",
			   				#"", "", "", ""],
			   			   #["" , "" ,
			   			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
			   			   "", "", "", ""],
			   			   #beta_3
			   			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			   			   #"", "", final_ests_point_all[3,1], "",
			   			   #"", "", "", ""],
			   			   #["" , "" ,
			   			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
			   			   "", "", "", ""],
			   			   # beta_4
			   			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			   			   #"", "", "", final_ests_point_all[4,1],
			   			   #"", "", "", ""],
			   			   #["" , "" ,
			   			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
			   			   "", "", "", ""],
			   			   #beta_5
			   			   ["Economies of scope", "", "", "", "", "", "", "", "", ""],
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
			   			   #"", "", "", "",
			   			   #final_ests_point_all[5,1], "", "", ""],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
			   			   #beta_6
			   			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
			   			   #"", "", "", "",
			   			   #"", final_ests_point_all[6,1], "", ""],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
			   			   #beta_7
			   			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
			   			   #"", "", "", "",
			   			   #"", "", final_ests_point_all[7,1], ""],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
			   			   #beta_8
			   			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
			   			   #"", "", "", "",
			   			   #"", "", "", final_ests_point_all[8,1]],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
			   			   #gamma merger cost (note the sign is reverse)
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   ["merger cost", "-\$\\gamma\$",
			   			   #-final_ests_point_all[1,2], -final_ests_point_all[2,2],
			   			   #-final_ests_point_all[3,2], -final_ests_point_all[4,2],
			   			   #-final_ests_point_all[5,2], -final_ests_point_all[6,2],
			   			   #-final_ests_point_all[7,2], -final_ests_point_all[8,2]],
			   			   #["" , "" ,
			   			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
			   			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
			   			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
			   			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
			   			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
			   			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
			   			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
			   			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
			   			   #delta subsidy sensitivity
			   			   ["subsidy sensitivity", L"\delta",
			   			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   Rule(),           # a nice \hline to make it ugly
			   			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all],
			   				#=["\\Correct Ineq" , "" ,
			    			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
			    			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
			    			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
			    			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
			    			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
			    			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
			    			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
			    			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
			   			   =#
			   			   ["\\% Inequalities" , "" ,
			   			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   			   Rule(:bottom)])




# check behavior
temp_subsidy_type = "to_buyer"
file_name_variable = "x1"
size_of_subsample_temp = 0
myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)

model_all_length = length(file_name_variable_list)
myests_point_all = zeros(model_all_length,
                         size(myests_point_scale_X_only)[1],
						 size(myests_point_scale_X_only)[2])
num_correct_ineq_all = zeros(model_all_length,
                         size(num_correct_ineq_scale_X_only)[1],
						 size(num_correct_ineq_scale_X_only)[2])
num_total_ineq_all = zeros(model_all_length,
                         size(num_total_ineq_scale_X_only)[1],
						 size(num_total_ineq_scale_X_only)[2])
accuracy_all = zeros(model_all_length,
                         size(accuracy_scale_X_only)[1])
final_ests_point_all = zeros(model_all_length,
                         size(final_ests_point_scale_X_only)[1])
for ii = 1:length(file_name_variable_list)
	temp_subsidy_type = "to_buyer"
	@show temp_file_name = file_name_variable_list[ii]
	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(num_correct_ineq_all[ii,:,1])[2], :],
	                                    digits=1)
end

CI_all_table = zeros(model_all_length,2,2)
for ii = 1:length(file_name_variable_list)
	temp_myests_point_all = myests_point_all[ii,num_correct_ineq_all[ii,:,1].==Int(findmax(num_correct_ineq_all[ii,:,1])[1]),:]
	CI_all_table[ii,:,:] = round.(hcat(
		Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
		Statistics.quantile(-temp_myests_point_all[:,2], [0.0,1.0])
		),digits=1)
end

total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy_main_firms_only.tex",
			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			  [Rule(:top),
			   ["","","", "", "", "", "", "", "", ""],
			   #["","",
			   #"Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   ["","",
			   "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]",
			    "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]"],
			   Rule(:mid),
			   ["Economies of scale", "", "", "", "", "", "", ""],
			   #beta_0
			   ["", "", "", "", "", "", "", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
			    "(S)", "(S)", "(S)", "(S)"],
			   #beta_1
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
			   #final_ests_point_all[1,1], "", "", "",
			   # "", "", "", ""],
			   #["" , "" ,
			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
			   "", "", "",
			   "", "", "", ""],
			   #beta_2
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
				#"", final_ests_point_all[2,1], "", "",
				#"", "", "", ""],
			   #["" , "" ,
			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
			   "", "", "", ""],
			   #beta_3
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			   #"", "", final_ests_point_all[3,1], "",
			   #"", "", "", ""],
			   #["" , "" ,
			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
			   "", "", "", ""],
			   # beta_4
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			   #"", "", "", final_ests_point_all[4,1],
			   #"", "", "", ""],
			   #["" , "" ,
			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
			   "", "", "", ""],
			   #beta_5
			   ["Economies of scope", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
			   #"", "", "", "",
			   #final_ests_point_all[5,1], "", "", ""],
			   #["" , "" ,
			   "", "", "", "",
			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
			   #beta_6
			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
			   #"", "", "", "",
			   #"", final_ests_point_all[6,1], "", ""],
			   #["" , "" ,
			   "", "", "", "",
			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
			   #beta_7
			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
			   #"", "", "", "",
			   #"", "", final_ests_point_all[7,1], ""],
			   #["" , "" ,
			   "", "", "", "",
			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
			   #beta_8
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
			   #"", "", "", "",
			   #"", "", "", final_ests_point_all[8,1]],
			   #["" , "" ,
			   "", "", "", "",
			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
			   #gamma merger cost (note the sign is reverse)
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", "", ""],
			   ["merger cost", "-\$\\gamma\$",
			   #-final_ests_point_all[1,2], -final_ests_point_all[2,2],
			   #-final_ests_point_all[3,2], -final_ests_point_all[4,2],
			   #-final_ests_point_all[5,2], -final_ests_point_all[6,2],
			   #-final_ests_point_all[7,2], -final_ests_point_all[8,2]],
			   #["" , "" ,
			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
			   #delta subsidy sensitivity
			   ["subsidy sensitivity", L"\delta",
			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
			   ["", "", "", "", "", "", "", "", "", ""],
			   Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all,
				total_ineq_all],
				#=["\\Correct Ineq" , "" ,
 			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
 			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
 			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
 			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
 			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
 			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
 			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
 			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
			   =#
			   ["\\% Inequalities" , "" ,
			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   Rule(:bottom)])


			   # check behavior
			   temp_subsidy_type = "shared"
			   file_name_variable = "x1"
			   size_of_subsample_temp = 0
			   myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
			   num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
			   num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
			   accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
			   final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
			   myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)

			   model_all_length = length(file_name_variable_list)
			   myests_point_all = zeros(model_all_length,
			                            size(myests_point_scale_X_only)[1],
			   						 size(myests_point_scale_X_only)[2])
			   num_correct_ineq_all = zeros(model_all_length,
			                            size(num_correct_ineq_scale_X_only)[1],
			   						 size(num_correct_ineq_scale_X_only)[2])
			   num_total_ineq_all = zeros(model_all_length,
			                            size(num_total_ineq_scale_X_only)[1],
			   						 size(num_total_ineq_scale_X_only)[2])
			   accuracy_all = zeros(model_all_length,
			                            size(accuracy_scale_X_only)[1])
			   final_ests_point_all = zeros(model_all_length,
			                            size(final_ests_point_scale_X_only)[1])
			   for ii = 1:length(file_name_variable_list)
			   	temp_subsidy_type = "shared"
			   	@show temp_file_name = file_name_variable_list[ii]
			   	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
			   	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
			   	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
			   	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
			   	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(num_correct_ineq_all[ii,:,1])[2], :],
			   	                                    digits=1)
			   end

			   CI_all_table = zeros(model_all_length,2,2)
			   for ii = 1:length(file_name_variable_list)
			   	temp_myests_point_all = myests_point_all[ii,num_correct_ineq_all[ii,:,1].==Int(findmax(num_correct_ineq_all[ii,:,1])[1]),:]
			   	CI_all_table[ii,:,:] = round.(hcat(
			   		Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
			   		Statistics.quantile(-temp_myests_point_all[:,2], [0.0,1.0])
			   		),digits=1)
			   end

			   total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
			   LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy_main_firms_only.tex",
			   			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			   			  [Rule(:top),
			   			   ["","","", "", "", "", "", "", "", ""],
			   			   #["","",
			   			   #"Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   			   ["","",
			   			   "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]",
			   			    "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]"],
			   			   Rule(:mid),
			   			   ["Economies of scale", "", "", "", "", "", "", ""],
			   			   #beta_0
			   			   ["", "", "", "", "", "", "", "", ""],
			   			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
			   			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
			   			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
			   			    "(S)", "(S)", "(S)", "(S)"],
			   			   #beta_1
			   			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
			   			   #final_ests_point_all[1,1], "", "", "",
			   			   # "", "", "", ""],
			   			   #["" , "" ,
			   			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
			   			   "", "", "",
			   			   "", "", "", ""],
			   			   #beta_2
			   			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
			   				#"", final_ests_point_all[2,1], "", "",
			   				#"", "", "", ""],
			   			   #["" , "" ,
			   			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
			   			   "", "", "", ""],
			   			   #beta_3
			   			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			   			   #"", "", final_ests_point_all[3,1], "",
			   			   #"", "", "", ""],
			   			   #["" , "" ,
			   			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
			   			   "", "", "", ""],
			   			   # beta_4
			   			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			   			   #"", "", "", final_ests_point_all[4,1],
			   			   #"", "", "", ""],
			   			   #["" , "" ,
			   			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
			   			   "", "", "", ""],
			   			   #beta_5
			   			   ["Economies of scope", "", "", "", "", "", "", "", "", ""],
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
			   			   #"", "", "", "",
			   			   #final_ests_point_all[5,1], "", "", ""],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
			   			   #beta_6
			   			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
			   			   #"", "", "", "",
			   			   #"", final_ests_point_all[6,1], "", ""],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
			   			   #beta_7
			   			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
			   			   #"", "", "", "",
			   			   #"", "", final_ests_point_all[7,1], ""],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
			   			   #beta_8
			   			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
			   			   #"", "", "", "",
			   			   #"", "", "", final_ests_point_all[8,1]],
			   			   #["" , "" ,
			   			   "", "", "", "",
			   			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
			   			   #gamma merger cost (note the sign is reverse)
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   ["merger cost", "-\$\\gamma\$",
			   			   #-final_ests_point_all[1,2], -final_ests_point_all[2,2],
			   			   #-final_ests_point_all[3,2], -final_ests_point_all[4,2],
			   			   #-final_ests_point_all[5,2], -final_ests_point_all[6,2],
			   			   #-final_ests_point_all[7,2], -final_ests_point_all[8,2]],
			   			   #["" , "" ,
			   			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
			   			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
			   			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
			   			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
			   			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
			   			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
			   			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
			   			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
			   			   #delta subsidy sensitivity
			   			   ["subsidy sensitivity", L"\delta",
			   			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
			   			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
			   			   ["", "", "", "", "", "", "", "", "", ""],
			   			   Rule(),           # a nice \hline to make it ugly
			   			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all,
			   				total_ineq_all],
			   				#=["\\Correct Ineq" , "" ,
			    			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
			    			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
			    			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
			    			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
			    			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
			    			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
			    			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
			    			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
			   			   =#
			   			   ["\\% Inequalities" , "" ,
			   			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   			   Rule(:bottom)])
