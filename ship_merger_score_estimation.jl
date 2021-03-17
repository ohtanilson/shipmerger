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
data = CSV.read("data_for_maximum_rank_estimation.csv")
rename!(data, :Column1 => :firm_id)

num_agents = size(data)[1]
names(data)
gdf = groupby(data, :group)
liner_sum = combine(gdf, [:liner] => sum)
special_sum = combine(gdf, [:special] => sum)
tanker_sum = combine(gdf, [:tanker] => sum)
tramper_sum = combine(gdf, [:tramper] => sum)

function extract_buyer_covariates_if_ii_is_buyer(ii, data)
	# pick individual own covariates
	buyer_X_scale = convert(Vector,data[ii,5:8])
	buyer_X_scope = convert(Vector,data[ii,11:14])
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
	deviater_X_scale = convert(Vector,subdata[kk,5:8])
	deviater_X_scope = convert(Vector,subdata[kk,11:14])
	res = vcat(deviater_X_scale, deviater_X_scope)
	return res
end

function extract_target_covariates_if_ii_is_buyer(ii, data)
	# pick up group names
	group_names = data[ii,4]
	# pick individual own covariates
	buyer_X_scale = convert(Vector,data[ii,5:8])
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

function extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(ii, kk, data)
	# pick up group names
	group_names = data[ii,4]
	# pick up deviating member
	subdata = @linq data |>
	    where(:group .== group_names)
	subdata.id[kk] == data.id[ii]
	deviater_X_scale = convert(Vector, subdata[kk,5:8])
	deviater_total_sum = subdata[kk,9]
	# pick individual own covariates
	buyer_X_scale = convert(Vector,data[ii,5:8])
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

function extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(ii, kk, data)
	group_names = data[ii,4]
	# pick individual own covariates
	buyer_X_scale = convert(Vector,data[ii,5:8])
	buyer_total_sum = data.total[ii]
	target_total_sum = data.group_total_tonnage[ii]
	# pick unmatched kk covariates
	unmatched_X_scale = convert(Vector,data[kk,5:8])
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
	merger_cost = num_of_firms_in_coalition/log(total_tonnage)
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

function gen_utility_est_without_subsidy(ii, buyer1_X::Vector, target1_X::Vector, data,
	                     beta, gamma, delta)
		#interaction_X_beta = beta.*log.(buyer1_X.*target1_X.+1)
		buyer1_X_total = sum(buyer1_X[1:4])
		target1_X_total = sum(target1_X[1:4])
		interaction_X_beta = 1.0*buyer1_X_total.*target1_X_total .+ beta.*buyer1_X.*target1_X
		#interaction_X_beta = beta.*buyer1_X.*target1_X
		merger_cost = gen_merger_cost_for_buyer_ii(ii, buyer1_X, target1_X, data)
		utility = sum(interaction_X_beta) - gamma*merger_cost
	return utility
end

function gen_unmatched_utility_est(buyer1_X::Vector, beta)
		#interaction_X_beta = beta.*log.(buyer1_X.+1)
		buyer1_X_total = sum(buyer1_X[1:4])
		#interaction_X_beta = 1.0*buyer1_X_total .+ beta.*buyer1_X
		interaction_X_beta = 1.0*buyer1_X_total .+ beta.*buyer1_X.*buyer1_X
		#interaction_X_beta = beta.*buyer1_X
		utility = sum(interaction_X_beta)
	return utility
end


function score_b_est_data(subsampled_id_list,
	                      data,
	                      theta::Vector{Float64},
						  subsidy_type)
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
			buyer2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data)
			# pick up target covariates
		    target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data)
			#First, I construct matching maximum score inequality for two acquiring firms without price data:
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target2_X, data, beta, gamma, delta, subsidy_type) # swapped
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target1_X, data, beta, gamma, delta, subsidy_type) # swapped
			#ineq[i, kk, hh] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			ineq[i, 1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
			# compare utility with and without subsidy
			# buyer 1
			payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
			ineq[i, 2] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
			# buyer 2
			payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
			ineq[i, 3] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
		#elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && IS_Sell(m,idx[2],target_bundle_id[2])
		elseif type_list_firm1[i] == "(1) main" && (type_list_firm2[i] != "(1) main" && type_list_firm2[i] != "unmatched")
			#println("iter $i = Case 2: firm 1 is a buyer and firm 2 is a seller.")
			#Second, I construct inequalities from an observed coalition:
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition of buyer's bundle
			subdata = @linq data |>
			    where(:group .== data[idx[1],4])
		    # construct swapped matches
			for kk = 1:size(subdata, 1)
				if subdata.id[kk] != data.id[idx[1]]
					target1_X_dropping_kk = extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(idx[1], kk, data)
					payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X_dropping_kk, data, beta, gamma, delta, subsidy_type)
					buyer2_X_unmatched_deviator_kk = extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(idx[1], kk, data)
					payoff_unobs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched_deviator_kk, beta)
					ineq[i,kk] = payoff_obs_match1 - payoff_unobs_match1 - payoff_unobs_match2
				else
					ineq[i,kk] = 0 # buyer and dropped firm are the same so ill-defined
					#println("buyer and dropped firm are the same so null valued.")
				end
			end
			# compare utility with and without subsidy
			# buyer 1
			payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
			ineq[i, size(subdata, 1)+1] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
		#elseif IS_Sell(m,idx[1],target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
	    elseif (type_list_firm1[i] != "(1) main" && type_list_firm1[i] != "unmatched") && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 3: firm 1 is a seller and firm 2 is a buyer.")
			#Second, I construct inequalities from an observed coalition:
			buyer2_X = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition of buyer's bundle
			subdata = @linq data |>
			    where(:group .== data[idx[2],4])
			# construct swapped matches
			for kk = 1:size(subdata, 1)
				if subdata.id[kk] != data.id[idx[2]]
					target2_X_dropping_kk = extract_target_covariates_if_ii_is_buyer_and_drop_member_id_kk(idx[2], kk, data)
					payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X_dropping_kk, data, beta, gamma, delta, subsidy_type)
					buyer1_X_unmatched_deviator_kk = extract_kk_buyer_covariates_if_if_ii_is_buyer_and_drop_member_id_kk(idx[2], kk, data)
					payoff_unobs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched_deviator_kk, beta)
					ineq[i,kk] = payoff_obs_match2 - payoff_unobs_match2 - payoff_unobs_match1
				else
				    ineq[i,kk] = 0 # buyer and dropped firm are the same so ill-defined
					#println("buyer and dropped firm are the same so null valued.")
				end
			end
			# compare utility with and without subsidy
			# buyer 2
			payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
			ineq[i, size(subdata, 1)+1] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
		#elseif IS_Buyer(m,idx[1],target_bundle_id[1]) && unmatched_vec == TB_toy(m.N, target_bundle_id[2])
	    elseif type_list_firm1[i] == "(1) main" && type_list_firm2[i] == "unmatched"
			#println("iter $i = Case 4: firm 1 is a buyer and firm 2 is unmatched.")
			#Third, I construct inequalities from an unmatched target:
			buyer1_X = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			target1_X = extract_target_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X, data, beta, gamma, delta, subsidy_type) # buyer
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_match2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta) # unmatched
			# construct swapped matches
			target1_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[1], idx[2], data)
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X, target1_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1
			# compare utility with and without subsidy
			# buyer 1
			payoff_unobs_match1_without_subsidy = gen_utility_est_without_subsidy(idx[1], buyer1_X, target1_X, data, beta, gamma, delta)
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_unmatch1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta)
			ineq[i,2] = payoff_obs_unmatch1 - payoff_unobs_match1_without_subsidy
		#elseif unmatched_vec == TB_toy(m.N, target_bundle_id[1]) && IS_Buyer(m,idx[2],target_bundle_id[2])
	    elseif type_list_firm1[i] == "unmatched" && type_list_firm2[i] == "(1) main"
			#println("iter $i = Case 5: firm 1 is unmatched and firm 2 is a buyer.")
			#Third, I construct inequalities from an unmatched target:
			buyer1_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[1], data)
			payoff_obs_match1 = gen_unmatched_utility_est(buyer1_X_unmatched, beta) # unmatched
			buyer2_X = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			target2_X = extract_target_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X, data, beta, gamma, delta, subsidy_type) # buyer
			# choose a firm out of coalition
			# construct swapped matches
			target2_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[2], idx[1], data)
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X, target2_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match2
			# buyer 2
			# compare utility with and without subsidy
			payoff_unobs_match2_without_subsidy = gen_utility_est_without_subsidy(idx[2], buyer2_X, target2_X, data, beta, gamma, delta)
			buyer2_X_unmatched = extract_buyer_covariates_if_ii_is_buyer(idx[2], data)
			payoff_obs_unmatch2 = gen_unmatched_utility_est(buyer2_X_unmatched, beta)
			ineq[i,2] = payoff_obs_unmatch2 - payoff_unobs_match2_without_subsidy
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
			target1_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[1], idx[2], data)
			payoff_unobs_match1 = gen_utility_est(idx[1], buyer1_X_unmatched, target1_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			target2_X_adding_unmatched_kk = extract_target_covariates_if_ii_is_buyer_and_adding_unmatched_kk(idx[2], idx[1], data)
			payoff_unobs_match2 = gen_utility_est(idx[2], buyer2_X_unmatched, target2_X_adding_unmatched_kk, data, beta, gamma, delta, subsidy_type) # swapped
			ineq[i,1,1] = payoff_obs_match1 + payoff_obs_match2 - payoff_unobs_match1 - payoff_unobs_match2
		else
			#println("iter $i = Case 7,8,9: no firms are buyers.")
			ineq[i, 1] = 0
		end
	end
	total_num_ineq = sum(ineq .!= 0)
	@show res = sum(ineq.>0)
	return res, total_num_ineq
end

#res,total_num_ineq = score_b_est_data(subsampled_id_list, data, theta)
function iterate_estimation(data,
	                        score_bthis::Function,
							subsidy_type
	                        ;param_dim=9,
	                        num_its=10,
	                        size_of_subsample=60,
							num_steps_DE=10)
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
		                                           size_of_subsample, replace = false)
		subsampled_id_list = vcat(main_firm_id, picked_not_main_firm_id)
		# Estimation
		m_res = BlackBoxOptim.bboptimize(theta -> score_bthis(subsampled_id_list,
		                                data, theta, subsidy_type)[1];
		                                SearchRange = (-100.0, 100.0),
										NumDimensions = param_dim,
										Method = :de_rand_1_bin,
										MaxSteps = num_steps_DE)
		#restore results
		println("score: ", score_bthis(subsampled_id_list, data, m_res.archive_output.best_candidate, subsidy_type)[2])
		~, num_correct_ineq[iter,1], num_total_ineq[iter,1] = score_bthis(subsampled_id_list,
		                                                                    data,
																			m_res.archive_output.best_candidate,
																			subsidy_type)
		println("the number of correct inequalities: ", num_correct_ineq[iter,1])
		println("the total number of inequalities: ", num_total_ineq[iter,1])
		myests[iter,:] = m_res.archive_output.best_candidate
	end
	return num_correct_ineq, num_total_ineq, myests
end

#=
function iterate_estimation_paralell(data,
	                        score_bthis::Function,
	                        ;param_dim=9,
	                        num_its=10,
	                        size_of_subsample=40,
							num_steps_DE=10)
	myests = zeros(num_its, param_dim)
	num_correct_ineq = zeros(num_its, 1)
	num_total_ineq = zeros(num_its, 1)
	data_main_firm = @linq data |>
		where(:type .== "(1) main")
	data_not_main_firm = @linq data |>
			where(:type .!= "(1) main")
	main_firm_id = data_main_firm.firm_id
	Threads.threading_run() do
	    i = threadid()
	    data[i] = test()
		picked_not_main_firm_id = StatsBase.sample(data_not_main_firm.firm_id,
		                                           size_of_subsample, replace = false)
		subsampled_id_list = vcat(main_firm_id, picked_not_main_firm_id)
		# Estimation
		m_res = BlackBoxOptim.bboptimize(theta -> score_bthis(subsampled_id_list, data, theta)[1];
		                                SearchRange = (-50.0, 50.0),
										NumDimensions = param_dim,
										Method = :de_rand_1_bin,
										MaxSteps = num_steps_DE)
		#restore results
		println("score: ", score_bthis(subsampled_id_list, data, m_res.archive_output.best_candidate)[2])
		~, num_correct_ineq[iter,1], num_total_ineq[iter,1] = score_bthis(subsampled_id_list,
		                                                                    data,
																			m_res.archive_output.best_candidate)
		println("the number of correct inequalities: ", num_correct_ineq[iter,1])
		println("the total number of inequalities: ", num_total_ineq[iter,1])
		myests[iter,:] = m_res.archive_output.best_candidate
	end
	for iter = 1:num_its
		Random.seed!(iter)
	end
	return num_correct_ineq, num_total_ineq, myests
end
=#

#point estimate
function score_bthis_full_X(subsampled_id_list, data, theta, subsidy_type)
	#target_theta = vcat(1,theta) # first parameter must be normalized to 1
    target_theta = vcat(theta) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data, target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end
function score_bthis_scale_X_only(subsampled_id_list, data, theta, subsidy_type)
	#target_theta = vcat(1, theta[1:3], zeros(4),theta[4:5]) # first parameter must be normalized to 1
	target_theta = vcat(theta[1:4], zeros(4),theta[5:6]) # first parameter must be normalized to 1
	score_res, total_num_ineq = score_b_est_data(subsampled_id_list, data, target_theta, subsidy_type)
	res = -1.0.*score_res .+ 100000.0 # need to be Float64 for bboptimize
	println("score:$(res) \n" )
	return res, score_res, total_num_ineq
end
function score_bthis_scale_and_scope_X_only(subsampled_id_list, data, theta, subsidy_type)
	#target_theta = vcat(1, zeros(3), theta) # first parameter must be normalized to 1
	target_theta = vcat(zeros(4), theta) # first parameter must be normalized to 1
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
                                                      param_dim = 10,
													  num_its=1, size_of_subsample=60,
													  num_steps_DE=2)
@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
                                                      score_bthis_scale_X_only,
													  temp_subsidy_type,
                                                      param_dim = 6,
													  num_its=1, size_of_subsample=60,
													  num_steps_DE=2)
@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
                                                      score_bthis_scale_and_scope_X_only,
													  temp_subsidy_type,
                                                      param_dim = 6,
													  num_its=1, size_of_subsample=60,
													  num_steps_DE=2)

#---------------------------#
# point estimate
#---------------------------#
num_steps_DE_temp = 100
num_its_temp = 1
size_of_fullsample = length(data.firm_id)-12
function point_estimate(subsidy_type;
	                    size_of_fullsample = 106,
	                    num_steps_DE_temp = 100,
	                    num_its_temp = 5)
	# model 1
    @time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                  score_bthis_scale_X_only,
													  subsidy_type,
                                                      param_dim = 6,
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
	                                                      param_dim = 6,
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
	                                                      param_dim = 10,
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
end

#----------------------------------------------------------#
# construct 95 percent Confidence interval via bootstrap
#----------------------------------------------------------#
function construct_CI(subsidy_type;
	                  num_its_bootstrap = 100,
	                  num_steps_DE_temp = 100,
	                  size_of_subsample_temp = 60)
    # model 1
	@time num_correct_ineq, num_total_ineq, myests = iterate_estimation(data,
	                                                  score_bthis_scale_X_only,
													  subsidy_type,
                                                      param_dim = 6,
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
	                                                      param_dim = 6,
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
	                                                      param_dim = 10,
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
you_want_run = "not run"
temp_subsidy_type = "shared"
if you_want_run == "run"
	@time point_estimate(temp_subsidy_type,
	                     num_steps_DE_temp = 100,
	                     #num_its_temp = 100,
						 num_its_temp = 2,
						 size_of_fullsample = 106)
	#@time point_estimate(num_steps_DE_temp = 100, num_its_temp = 1)
	#246.815583 seconds (1.28 G allocations: 124.312 GiB, 5.89% gc time)
	#@time point_estimate(num_steps_DE_temp = 100,  num_its_temp = 10)
	#3812.250026 seconds (20.37 G allocations: 2.009 TiB, 6.69% gc time)
	@time construct_CI(temp_subsidy_type,
	                   num_its_bootstrap = 2,
	                   num_steps_DE_temp = 100,
					   size_of_subsample_temp = 60)
	#34558.137605 seconds (162.34 G allocations: 15.437 TiB, 5.96% gc time)
end

you_want_run = "not run"
temp_subsidy_type = "to_buyer"
if you_want_run == "run"
	@time point_estimate(temp_subsidy_type,
	                     num_steps_DE_temp = 100,
	                     #num_its_temp = 100,
						 num_its_temp = 100,
						 size_of_fullsample = 106)
	#@time point_estimate(num_steps_DE_temp = 100, num_its_temp = 1)
	#246.815583 seconds (1.28 G allocations: 124.312 GiB, 5.89% gc time)
	#@time point_estimate(num_steps_DE_temp = 100,  num_its_temp = 10)
	#3812.250026 seconds (20.37 G allocations: 2.009 TiB, 6.69% gc time)
	@time construct_CI(temp_subsidy_type,
	                   num_its_bootstrap = 100,
	                   num_steps_DE_temp = 100,
					   size_of_subsample_temp = 60)
	#34558.137605 seconds (162.34 G allocations: 15.437 TiB, 5.96% gc time)
end
#--------------------#
# read output (shared subsidy)
#--------------------#
# model 1
temp_subsidy_type = "shared"
size_of_subsample_temp = 60
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

# CI
#beta = theta[1:8]
#gamma = theta[9] # coefficient on merger cost
#delta = theta[10] # coefficient on subsidy indicator
myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
myests_CI_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
myests_CI_full_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
CI_scale_X_only_table = round.(hcat(
    Statistics.quantile(myests_CI_scale_X_only[:,1], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,2], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,3], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,4], [0.025,0.975]),
	Statistics.quantile(-myests_CI_scale_X_only[:,5], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,6], [0.025,0.975])
	),digits=2)
CI_scope_X_only_table = round.(hcat(
    Statistics.quantile(myests_CI_scope_X_only[:,1], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,2], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,3], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,4], [0.025,0.975]),
	Statistics.quantile(-myests_CI_scope_X_only[:,5], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,6], [0.025,0.975])
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
	Statistics.quantile(-myests_CI_full_X_only[:,9], [0.025,0.975]),
	Statistics.quantile(myests_CI_full_X_only[:,10], [0.025,0.975])
	),digits=2)

LaTeXTabulars.latex_tabular("ship_merger/figuretable/score_results_temp_$(temp_subsidy_type)_subsidy.tex",
              Tabular("@{\\extracolsep{5pt}}lcccc"),
              [Rule(:top),
			   ["","","", "Value Function", ""],
               ["","","Point Estimate", "Point Estimate", "Point Estimate"],
			   ["","","[95\\% CI Set Identified]", "[95\\% CI Set Identified]", "[95\\% CI Set Identified]"],
			   Rule(:mid),

			   ["Measure of economies of scale", "", "", "", ""],
			   #beta_0
			   ["", "", "", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0", "+1", "+1", "+1"],
			   ["" , "" , "(Superconsistent)", "(Superconsistent)", "(Superconsistent)"],
			   #beta_1
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
			    final_ests_point_scale_X_only[1], "", final_ests_point_full_X_only[1]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,1]), $(CI_scale_X_only_table[2,1])]",
			   "",
			   "[$(CI_full_X_only_table[1,1]), $(CI_full_X_only_table[2,1])]"],
			   #beta_2
               ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
			    final_ests_point_scale_X_only[2], "", final_ests_point_full_X_only[2]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,2]), $(CI_scale_X_only_table[2,2])]",
			   "",
			   "[$(CI_full_X_only_table[1,2]), $(CI_full_X_only_table[2,2])]"],
			   #beta_3
               ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			    final_ests_point_scale_X_only[3], "", final_ests_point_full_X_only[3]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,3]), $(CI_scale_X_only_table[2,3])]",
			   "",
			   "[$(CI_full_X_only_table[1,3]), $(CI_full_X_only_table[2,3])]"],
			   # beta_4
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			    final_ests_point_scale_X_only[4], "", final_ests_point_full_X_only[4]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,4]), $(CI_scale_X_only_table[2,4])]",
			   "",
			   "[$(CI_full_X_only_table[1,4]), $(CI_full_X_only_table[2,4])]"],
			   #beta_5
			   ["Measure of economies of scope", "", "", "", ""],
			   ["", "", "", "", ""],
               ["share of liner\$_{b}\$ \$\\times\$ share of liner\$_{t}\$", L"\beta_5",
			    "", final_ests_point_scope_X_only[1], final_ests_point_full_X_only[5]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,1]), $(CI_scope_X_only_table[2,1])]",
			   "[$(CI_full_X_only_table[1,5]), $(CI_full_X_only_table[2,5])]"],
			   #beta_6
			   ["share of tramper\$_{b}\$ \$\\times\$ share of tramper\$_{t}\$", L"\beta_6",
			    "", final_ests_point_scope_X_only[2], final_ests_point_full_X_only[6]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,2]), $(CI_scope_X_only_table[2,2])]",
			   "[$(CI_full_X_only_table[1,6]), $(CI_full_X_only_table[2,6])]"],
			   #beta_7
               ["share of special\$_{b}\$ \$\\times\$ share of special\$_{t}\$", L"\beta_7",
			    "", final_ests_point_scope_X_only[3], final_ests_point_full_X_only[7]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,3]), $(CI_scope_X_only_table[2,3])]",
			   "[$(CI_full_X_only_table[1,7]), $(CI_full_X_only_table[2,7])]"],
			   #beta_8
               ["share of tanker\$_{b}\$ \$\\times\$ share of tanker\$_{t}\$", L"\beta_8",
			   "", final_ests_point_scope_X_only[4], final_ests_point_full_X_only[8]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,4]), $(CI_scope_X_only_table[2,4])]",
			   "[$(CI_full_X_only_table[1,8]), $(CI_full_X_only_table[2,8])]"],
			   #gamma merger cost (note the sign is reverse)
			   ["", "", "", "", ""],
			   ["", "", "", "", ""],
               ["merger cost", L"\gamma",
			   -final_ests_point_scale_X_only[5], -final_ests_point_scope_X_only[5], -final_ests_point_full_X_only[9]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,4]), $(CI_scale_X_only_table[2,4])]",
			   "[$(CI_scope_X_only_table[1,5]), $(CI_scope_X_only_table[2,5])]",
			   "[$(CI_full_X_only_table[1,9]), $(CI_full_X_only_table[2,9])]"],
			   #delta subsidy sensitivity
               ["subsidy sensitivity (\$s^\{\\text\{shared\}\}\$)", L"\delta",
			   final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,5]), $(CI_scale_X_only_table[2,5])]",
			   "[$(CI_scope_X_only_table[1,6]), $(CI_scope_X_only_table[2,6])]",
			   "[$(CI_full_X_only_table[1,10]), $(CI_full_X_only_table[2,10])]"],
               Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities in Point Estimate" , "" ,
			    Int64(num_total_ineq_scale_X_only[findmax(accuracy_scale_X_only)[2]]),
				Int64(num_total_ineq_scope_X_only[findmax(accuracy_scope_X_only)[2]]),
				Int64(num_total_ineq_full_X_only[findmax(accuracy_full_X_only)[2]])],
			   ["\\% Inequalities" , "" ,
			   round(accuracy_scale_X_only[findmax(accuracy_scale_X_only)[2]],digits=3),
			   round(accuracy_scope_X_only[findmax(accuracy_scope_X_only)[2]],digits=3),
			   round(accuracy_full_X_only[findmax(accuracy_full_X_only)[2]],digits=3)],
               Rule(:bottom)])

#--------------------#
# read output (to_buyer subsidy)
#--------------------#
# model 1
temp_subsidy_type = "to_buyer"
size_of_subsample_temp = 60
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

# CI
#beta = theta[1:8]
#gamma = theta[9] # coefficient on merger cost
#delta = theta[10] # coefficient on subsidy indicator
myests_CI_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
myests_CI_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
myests_CI_full_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
CI_scale_X_only_table = round.(hcat(
    Statistics.quantile(myests_CI_scale_X_only[:,1], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,2], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,3], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,4], [0.025,0.975]),
	Statistics.quantile(-myests_CI_scale_X_only[:,5], [0.025,0.975]),
	Statistics.quantile(myests_CI_scale_X_only[:,6], [0.025,0.975])
	),digits=2)
CI_scope_X_only_table = round.(hcat(
    Statistics.quantile(myests_CI_scope_X_only[:,1], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,2], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,3], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,4], [0.025,0.975]),
	Statistics.quantile(-myests_CI_scope_X_only[:,5], [0.025,0.975]),
	Statistics.quantile(myests_CI_scope_X_only[:,6], [0.025,0.975])
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
	Statistics.quantile(-myests_CI_full_X_only[:,9], [0.025,0.975]),
	Statistics.quantile(myests_CI_full_X_only[:,10], [0.025,0.975])
	),digits=2)

LaTeXTabulars.latex_tabular("ship_merger/figuretable/score_results_temp_$(temp_subsidy_type)_subsidy.tex",
              Tabular("@{\\extracolsep{5pt}}lcccc"),
              [Rule(:top),
			   ["","","", "Value Function", ""],
               ["","","Point Estimate", "Point Estimate", "Point Estimate"],
			   ["","","[95\\% CI Set Identified]", "[95\\% CI Set Identified]", "[95\\% CI Set Identified]"],
			   Rule(:mid),

			   ["Measure of economies of scale", "", "", "", ""],
			   #beta_0
			   ["", "", "", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0", "+1", "+1", "+1"],
			   ["" , "" , "(Superconsistent)", "(Superconsistent)", "(Superconsistent)"],
			   #beta_1
			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
			    final_ests_point_scale_X_only[1], "", final_ests_point_full_X_only[1]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,1]), $(CI_scale_X_only_table[2,1])]",
			   "",
			   "[$(CI_full_X_only_table[1,1]), $(CI_full_X_only_table[2,1])]"],
			   #beta_2
               ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
			    final_ests_point_scale_X_only[2], "", final_ests_point_full_X_only[2]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,2]), $(CI_scale_X_only_table[2,2])]",
			   "",
			   "[$(CI_full_X_only_table[1,2]), $(CI_full_X_only_table[2,2])]"],
			   #beta_3
               ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
			    final_ests_point_scale_X_only[3], "", final_ests_point_full_X_only[3]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,3]), $(CI_scale_X_only_table[2,3])]",
			   "",
			   "[$(CI_full_X_only_table[1,3]), $(CI_full_X_only_table[2,3])]"],
			   # beta_4
			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
			    final_ests_point_scale_X_only[4], "", final_ests_point_full_X_only[4]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,4]), $(CI_scale_X_only_table[2,4])]",
			   "",
			   "[$(CI_full_X_only_table[1,4]), $(CI_full_X_only_table[2,4])]"],
			   #beta_5
			   ["Measure of economies of scope", "", "", "", ""],
			   ["", "", "", "", ""],
               ["share of liner\$_{b}\$ \$\\times\$ share of liner\$_{t}\$", L"\beta_5",
			    "", final_ests_point_scope_X_only[1], final_ests_point_full_X_only[5]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,1]), $(CI_scope_X_only_table[2,1])]",
			   "[$(CI_full_X_only_table[1,5]), $(CI_full_X_only_table[2,5])]"],
			   #beta_6
			   ["share of tramper\$_{b}\$ \$\\times\$ share of tramper\$_{t}\$", L"\beta_6",
			    "", final_ests_point_scope_X_only[2], final_ests_point_full_X_only[6]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,2]), $(CI_scope_X_only_table[2,2])]",
			   "[$(CI_full_X_only_table[1,6]), $(CI_full_X_only_table[2,6])]"],
			   #beta_7
               ["share of special\$_{b}\$ \$\\times\$ share of special\$_{t}\$", L"\beta_7",
			    "", final_ests_point_scope_X_only[3], final_ests_point_full_X_only[7]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,3]), $(CI_scope_X_only_table[2,3])]",
			   "[$(CI_full_X_only_table[1,7]), $(CI_full_X_only_table[2,7])]"],
			   #beta_8
               ["share of tanker\$_{b}\$ \$\\times\$ share of tanker\$_{t}\$", L"\beta_8",
			   "", final_ests_point_scope_X_only[4], final_ests_point_full_X_only[8]],
			   ["" , "" ,
			   "",
			   "[$(CI_scope_X_only_table[1,4]), $(CI_scope_X_only_table[2,4])]",
			   "[$(CI_full_X_only_table[1,8]), $(CI_full_X_only_table[2,8])]"],
			   #gamma merger cost (note the sign is reverse)
			   ["", "", "", "", ""],
			   ["", "", "", "", ""],
               ["merger cost", L"\gamma",
			   -final_ests_point_scale_X_only[5], -final_ests_point_scope_X_only[5], -final_ests_point_full_X_only[9]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,4]), $(CI_scale_X_only_table[2,4])]",
			   "[$(CI_scope_X_only_table[1,5]), $(CI_scope_X_only_table[2,5])]",
			   "[$(CI_full_X_only_table[1,9]), $(CI_full_X_only_table[2,9])]"],
			   #delta subsidy sensitivity
               ["subsidy sensitivity (\$s_{i,J}^{\\text{to buyer}}\$)", L"\delta",
			   final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   ["" , "" ,
			   "[$(CI_scale_X_only_table[1,5]), $(CI_scale_X_only_table[2,5])]",
			   "[$(CI_scope_X_only_table[1,6]), $(CI_scope_X_only_table[2,6])]",
			   "[$(CI_full_X_only_table[1,10]), $(CI_full_X_only_table[2,10])]"],
               Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities in Point Estimate" , "" ,
			    Int64(num_total_ineq_scale_X_only[findmax(accuracy_scale_X_only)[2]]),
				Int64(num_total_ineq_scope_X_only[findmax(accuracy_scope_X_only)[2]]),
				Int64(num_total_ineq_full_X_only[findmax(accuracy_full_X_only)[2]])],
			   ["\\% Inequalities" , "" ,
			   round(accuracy_scale_X_only[findmax(accuracy_scale_X_only)[2]],digits=3),
			   round(accuracy_scope_X_only[findmax(accuracy_scope_X_only)[2]],digits=3),
			   round(accuracy_full_X_only[findmax(accuracy_full_X_only)[2]],digits=3)],
               Rule(:bottom)])

#-----------------------------------------------------#
# find point-estimate LB of model 2 (shared subsidy)
#-----------------------------------------------------#
temp_subsidy_type = "shared"
myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
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
                xlabel = "β₅",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₅",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[1]],
             label = "β₅ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[1]],
             label = "β₅ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[1]],
             label = "β₅ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_5_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,2], bins=40,
                ylabel = "count",
                xlabel = "β₆",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₆",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[2]],
             label = "β₆ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[2]],
             label = "β₆ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[2]],
             label = "β₆ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_6_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,3], bins=40,
                ylabel = "count",
                xlabel = "β₇",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₇",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[3]],
             label = "β₇ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[3]],
             label = "β₇ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[3]],
             label = "β₇ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_7_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,4], bins=40,
                ylabel = "count",
                xlabel = "β₈",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₈",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[4]],
             label = "β₈ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[4]],
             label = "β₈ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[4]],
             label = "β₈ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_8_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(-myests_point_scope_X_only[:,5], bins=40,
                ylabel = "count",
                xlabel = "γ (merger cost)",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated γ",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([-final_ests_point_scope_X_only[5]],
             label = "γ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([-second_best_ests_point_scope_X_only[5]],
             label = "γ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([-third_best_ests_point_scope_X_only[5]],
             label = "γ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_delta_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,6], bins=40,
                ylabel = "count",
                xlabel = "δ (subsidy sensitivity)",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated δ",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[6]],
             label = "δ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[6]],
             label = "δ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[6]],
             label = "δ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_gamma_scope_X_only_$(temp_subsidy_type)_subsidy")

#--------------------------------------#
# find LB of model 2 (to buyer subsidy)
#--------------------------------------#
temp_subsidy_type = "to_buyer"

myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_106_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
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
                xlabel = "β₅",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₅",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[1]],
             label = "β₅ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[1]],
             label = "β₅ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[1]],
             label = "β₅ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_5_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,2], bins=40,
                ylabel = "count",
                xlabel = "β₆",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₆",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[2]],
             label = "β₆ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[2]],
             label = "β₆ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[2]],
             label = "β₆ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_6_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,3], bins=40,
                ylabel = "count",
                xlabel = "β₇",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₇",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[3]],
             label = "β₇ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[3]],
             label = "β₇ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[3]],
             label = "β₇ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_7_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,4], bins=40,
                ylabel = "count",
                xlabel = "β₈",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated β₈",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[4]],
             label = "β₈ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[4]],
             label = "β₈ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[4]],
             label = "β₈ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_beta_8_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(-myests_point_scope_X_only[:,5], bins=40,
                ylabel = "count",
                xlabel = "γ (merger cost)",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated γ",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([-final_ests_point_scope_X_only[5]],
             label = "γ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([-second_best_ests_point_scope_X_only[5]],
             label = "γ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([-third_best_ests_point_scope_X_only[5]],
             label = "γ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_delta_scope_X_only_$(temp_subsidy_type)_subsidy")
Plots.histogram(myests_point_scope_X_only[:,6], bins=40,
                ylabel = "count",
                xlabel = "δ (subsidy sensitivity)",
				title = "Subsidy type ($temp_subsidy_type)",
                label = "point-estimated δ",
				xlim =[-100,100],
				alpha = 0.3)
Plots.vline!([final_ests_point_scope_X_only[6]],
             label = "δ (highest score:$(best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "black")
Plots.vline!([second_best_ests_point_scope_X_only[6]],
             label = "δ (second best score:$(second_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "red")
Plots.vline!([third_best_ests_point_scope_X_only[6]],
             label = "δ (third best score:$(third_best_accuracy_scope_X_only))",
			 linestyle = :dash,
             color = "blue")
savefig("julia_merger_figure/histogram_estimated_parameters_gamma_scope_X_only_$(temp_subsidy_type)_subsidy")
