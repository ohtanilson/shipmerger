function score_b_est_data_only_HHI(subsampled_id_list,
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
	# theta[9] assigns coefficient on HHI
	beta = theta[1:9]#theta[1:8]
	gamma = theta[10]#theta[9] # coefficient on merger cost
	delta = theta[11]#theta[10] # coefficient on subsidy indicator
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
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
			buyer2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum, HHI_included = true)
			# pick up target covariates
		    target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data, info_sum = temp_info_sum, HHI_included = true)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum, HHI_included = true)
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
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, 2] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
				ineq[i, 3] = payoff_obs_match1 - payoff_obs_unmatch1 # with subsidy
				# buyer 2
				payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, 4] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
				ineq[i, 5] = payoff_obs_match2 - payoff_obs_unmatch2 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, 2] = payoff_obs_match1 - payoff_obs_unmatch1
				ineq[i, 3] = payoff_obs_match2 - payoff_obs_unmatch2
			end
		#elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && IS_Sell(m,idx[2],target_bundle_id[2])
		elseif type_list_firm1[i] == "(1) main" && (type_list_firm2[i] != "(1) main" && type_list_firm2[i] != "unmatched")
			#println("iter $i = Case 2: firm 1 is a buyer and firm 2 is a seller.")
			#Second, I construct inequalities from an observed coalition:
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
			target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition of buyer's bundle
			subdata = @linq data |>
			    where(:group .== data[idx[1],4])
		    # construct swapped matches
			for kk = 1:size(subdata, 1)
				if subdata.id[kk] != data.id[idx[1]]
					target1_X_dropping_kk = extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(idx[1], kk, data, info_sum = temp_info_sum, HHI_included = true)
					payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X_dropping_kk, data, beta, gamma, delta, subsidy_type)
					buyer2_X_unmatched_deviator_kk = extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(idx[1], kk, data, HHI_included = true)
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
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, size(subdata, 1)+1] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
				ineq[i, size(subdata, 1)+2] = payoff_obs_match1 - payoff_obs_unmatch1 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, size(subdata, 1)+3] = payoff_obs_match1 - payoff_obs_unmatch1
			end
		#elseif IS_Sell(m,idx[1],target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
	    elseif (type_list_firm1[i] != "(1) main" && type_list_firm1[i] != "unmatched") && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 3: firm 1 is a seller and firm 2 is a buyer.")
			#Second, I construct inequalities from an observed coalition:
			buyer2_X = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition of buyer's bundle
			subdata = @linq data |>
			    where(:group .== data[idx[2],4])
			# construct swapped matches
			for kk = 1:size(subdata, 1)
				if subdata.id[kk] != data.id[idx[2]]
					target2_X_dropping_kk = extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(idx[2], kk, data, info_sum = temp_info_sum, HHI_included = true)
					payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X_dropping_kk, data, beta, gamma, delta, subsidy_type)
					buyer1_X_unmatched_deviator_kk = extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(idx[2], kk, data, HHI_included = true)
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
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, size(subdata, 1)+1] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
				ineq[i, size(subdata, 1)+2] = payoff_obs_match2 - payoff_obs_unmatch2 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i, size(subdata, 1)+3] = payoff_obs_match2 - payoff_obs_unmatch2
			end
		#elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
	    elseif type_list_firm1[i] == "(1) main" && type_list_firm2[i] == "unmatched"
			#println("iter $i = Case 4: firm 1 is a buyer and firm 2 is unmatched.")
			#Third, I construct inequalities from an unmatched target:
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
			target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
			payoff_obs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta) # unmatched
			# construct swapped matches
			target1_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[1], idx[2], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i, 1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				# buyer 1
				payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, 2] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
				ineq[i, 3] = payoff_obs_match1 - payoff_obs_unmatch1 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
				payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
				ineq[i, 4] = payoff_obs_match1 - payoff_obs_unmatch1
			end
		#elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
	    elseif type_list_firm1[i] == "unmatched" && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 5: firm 1 is unmatched and firm 2 is a buyer.")
			#Third, I construct inequalities from an unmatched target:
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
			payoff_obs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta) # unmatched
			buyer2_X = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition
			# construct swapped matches
			target2_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[2], idx[1], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match2
			# buyer 2
			if compare_with_and_without_subsidy == "yes"
				# compare utility with and without subsidy
				payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i,2] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
				ineq[i,3] = payoff_obs_match2 - payoff_obs_unmatch2 # with subsidy
			end
			if compare_matched_and_unmatched == "yes"
				# compare matched utility with unmatched utility
				buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
				payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
				ineq[i,4] = payoff_obs_match2 - payoff_obs_unmatch2
			end
		#elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
	    elseif type_list_firm1[i] == "unmatched" && type_list_firm2[i] == "unmatched"
			#println("iter $i = Case 6: both picked firms are unmatched.")
			# Fourth, I construct inequalities from IR conditions without subsidy
			# Here, deviation comes from merging unmatched pairwise firm
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data, HHI_included = true)
			payoff_obs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta) # unmatched
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data, HHI_included = true)
			payoff_obs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta) # unmatched
			# construct swapped matches
			target1_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[1], idx[2], data, info_sum = temp_info_sum, HHI_included = true)
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X_unmatched, target1_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			target2_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[2], idx[1], data, info_sum = temp_info_sum, HHI_included = true)
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

function score_bthis_only_HHI(subsampled_id_list, data, theta, subsidy_type;
	                        calibrated_delta = 1)
	#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
	target_theta = vcat(zeros(4),
	                    theta[1],# HHI coefficient
						zeros(4),
						theta[2],
	                    calibrated_delta) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data_only_HHI(subsampled_id_list, data,
	                                                      target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end


function iterate_estimation_only_HHI(data,
	                        score_bthis::Function,
							subsidy_type
	                        ;param_dim=9,
	                        num_its=10,
	                        size_of_subsample=60,
							num_steps_DE=10,
							temp_calibrated_delta = 1,
							info_sum = temp_info_sum,
							sampling_method = "bootstrap")
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
		if sampling_method == "bootstrap"
		    picked_not_main_firm_id = StatsBase.sample(data_not_main_firm.firm_id,
		                                           size_of_subsample,
												   replace = true#false
												   )
		elseif sampling_method == "subsampling"
			picked_not_main_firm_id = StatsBase.sample(data_not_main_firm.firm_id,
		                                           size_of_subsample,
												   replace = false
												   )
		end
		subsampled_id_list = vcat(main_firm_id, picked_not_main_firm_id)
		# Estimation
		temp_search_domain = [(-200.0,200.0)]
		m_res = BlackBoxOptim.bboptimize(theta -> score_bthis(subsampled_id_list, data, theta, subsidy_type,
										              calibrated_delta = temp_calibrated_delta)[1];
		                                SearchRange = vcat(temp_search_domain[1],(-100.0, 2000.0)),
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

#---------------------------#
# point estimate
#---------------------------#
function point_estimate_only_HHI(subsidy_type;
	                    size_of_fullsample = 106,
	                    num_steps_DE_temp = 50,
	                    num_its_temp = 100,
						temp_temp_calibrated_delta = 1,
						temp_sampling_method = "subsampling")
	# model 1
    @time num_correct_ineq, num_total_ineq, myests = iterate_estimation_only_HHI(data,
	                                                  score_bthis_only_HHI,
													  subsidy_type,
													  temp_calibrated_delta = temp_temp_calibrated_delta,
                                                      param_dim = 2,
													  num_its = num_its_temp,
													  size_of_subsample = size_of_fullsample,
													  num_steps_DE = num_steps_DE_temp,
													  sampling_method = temp_sampling_method)
	open("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_HHI_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_HHI_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(subsidy_type)_subsidy_HHI_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end

end

#----------------------------------------------------------#
# construct 95 percent Confidence interval via bootstrap
#----------------------------------------------------------#
function construct_CI_only_HHI(subsidy_type;
	                  num_its_bootstrap = 100,
	                  num_steps_DE_temp = 50,
	                  size_of_subsample_temp = 60,
					  temp_temp_calibrated_delta = 1,
					  temp_sampling_method = "subsampling")
    # model 1
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation_only_HHI(data,
	                                                  score_bthis_only_HHI,
													  subsidy_type,
													  temp_calibrated_delta = temp_temp_calibrated_delta,
                                                      param_dim = 2,
													  num_its = num_its_bootstrap,
													  size_of_subsample = size_of_fullsample,
													  num_steps_DE = num_steps_DE_temp,
													  sampling_method = temp_sampling_method)
	open("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_HHI_only.txt", "w") do io
		DelimitedFiles.writedlm(io, myests,",")
	end
	open("julia_merger_result/CI_num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_HHI_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_total_ineq,",")
	end
	open("julia_merger_result/CI_num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(subsidy_type)_subsidy_HHI_only.txt", "w") do io
		DelimitedFiles.writedlm(io, num_correct_ineq,",")
	end
end
