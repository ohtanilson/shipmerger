model_id_iter = 1 # middle scenario
unique_matches_list = zeros(
                      length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list),
					  20, # preassign upper bound of unique number of matching outcomes
					  m.N,
					  m.PN)
println("find allocation indices different from benchmark")
different_index = zeros(length(threshold_tonnage_list),
                        length(subsidy_amount_list),
						iter_end,
						iter_end)
@time for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
	for kk_benchmark = 1:iter_end
		for kk = 1:iter_end
			# compare benchmark allocation (iter=1)
		    unique_index = sum(abs.(matches_list[model_id_iter, nn, mm, kk_benchmark, :, :].-matches_list[model_id_iter, nn, mm, kk, :, :]))
			if unique_index != 0
				# @show nn,mm,kk
				different_index[nn,mm,kk_benchmark,kk] = 1
				# unique_matches_list[model_id_iter, nn, mm, kk, :, :] = matches_list[model_id_iter, nn, mm, kk, :, :]
			end
		end
	end
end
most_frequently_observed_allcation_index = zeros(length(threshold_tonnage_list),
                                                 length(subsidy_amount_list))
@time for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
	@show nn,mm
    iter_min = findmin(sum(different_index[nn,mm,:,:],dims = 2))[2][1]
	@show most_frequently_observed_allcation_index[nn,mm] = iter_min
end
# assign most frequently observed allocations
target_matches_list = zeros(
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list),
					  m.N,
					  m.PN)
@time for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
    target_matches_list[nn,mm,:,:] = matches_list[model_id_iter,nn,mm,
	                                          Int.(most_frequently_observed_allcation_index[nn,mm]),
											  :,:]
end
# Case 1: actual amount with different thresholds
println("*****\nCase 1: actual amount with different thresholds\n******")
temp_matches = target_matches_list[:,6,:,:] # actual amount scenario
for temp_threshold = 1:length(threshold_tonnage_list)
	println("*****threshold $(threshold_tonnage_list[temp_threshold])******************")
	for ii = 1:m.N
		m.Bundle[temp_matches[temp_threshold,ii,:].>0]
		if m.Bundle[temp_matches[temp_threshold,ii,:].>0] != "000000000000"
			maxbundle = findmax(temp_matches[temp_threshold,ii,:])[2]
		    println("firm $ii buys $(m.Bundle[maxbundle])")
		else
			println("firm $ii is unmatched")
		end
	end
end
#-----------------------#
# interpretable results
#-----------------------#
println("
threshold_tonnage_list[temp_threshold] = 1.0
1 buys 4
2 buys 9
5 buys 10
6 buys 3
11 buys 7
12 buys 8
")

println("
threshold_tonnage_list[temp_threshold] = 2.0
1 buys 4
5 buys 7,10
6 buys 2,3,11
12 buys 2(duplicated),8,9
2 is duplicated because the allocation is 0.5 equally
")

println("
threshold_tonnage_list[temp_threshold] = 3.0
1 buys 7,10
5 buys 2,4,9,12
6 buys 2(duplicated),3,8,10
11 sells to whom?
2 is duplicated because the allocation is 0.5 equally
")

println("
threshold_tonnage_list[temp_threshold] = 4.0
1 buys 3,4,7,10
6 buys 2,3(duplicated),5,8,9,12
11 unmatched
3 is duplicated because the allocation is 0.5 equally
")

println("
threshold_tonnage_list[temp_threshold] = 5.0
1 buys 2,3,4,6,7,8,9,10
11 unmatched
12 unmatched
")

println("
threshold_tonnage_list[temp_threshold] = 7.5
1 buys 2,3,4,5,6,7,8,9,10,11,12
")

# Case 2: different amount with an actual threshold
println("*****\nCase 2: different amount with an actual threshold\n******************")
temp_matches = target_matches_list[1,:,:,:] # actual threshold scenario
for temp_amount = 1:length(subsidy_amount_list)
	println("******amount = $(subsidy_amount_list[temp_amount])*****************")
	subsidy_amount_list[temp_amount]
	for ii = 1:m.N
		m.Bundle[temp_matches[temp_amount,ii,:].>0]
		if m.Bundle[temp_matches[temp_amount,ii,:].>0] != "000000000000"
			maxbundle = findmax(temp_matches[temp_amount,ii,:])[2]
		    println("firm $ii buys $(m.Bundle[maxbundle])")
		else
			println("firm $ii is unmatched")
		end
	end
end
#-----------------------#
# interpretable results
#-----------------------#
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
2 buys 9
4 unmatched
5 buys 10
6 buys 3
11 buys 7
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 0.75
1 buys 4
2 buys 9
5 buys 10
6 buys 3
11 buys 7
12 buys 8
")
println("
subsidy_amount_list[temp_amount] = 1.0
1 buys 4
2 buys 9
5 buys 10
6 buys 3
11 buys 7
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 2.0
1 buys 4
2 buys 9
5 buys 10
6 buys 7
11 buys 3
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 3.0
1 buys 4
2 buys 9
5 buys 10
6 buys 7
11 buys 3
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 4.0
1 buys 4
2 buys 9
5 buys 10
6 buys 7
11 buys 3
12 buys 8
")

println("
subsidy_amount_list[temp_amount] = 5.0
1 buys 4
2 buys 9
5 buys 10
6 buys 7
11 buys 3
12 buys 8
")
