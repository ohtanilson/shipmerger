# this script must be located in ship_merger_score_estimation.jl
# check behavior
#temp_subsidy_type = "shared"
file_name_variable = "x1"
size_of_subsample_temp = 0 # main firms only
myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)

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
	#temp_subsidy_type = "shared"
	@show temp_file_name = file_name_variable_list[ii]
	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
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
LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy_main_firms_only_delta_$(temp_calibrated_delta_list[1]).tex",
			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
			  [Rule(:top),
			   ["","","", "", "", "", "", "", "", ""],
			   #["","",
			   #"Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   ["","",
			   "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]",
			    "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]"],
			   Rule(:mid),
			   ["Scale variables", "", "", "", "", "", "", ""],
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
			   ["Share variables", "", "", "", "", "", "", "", "", ""],
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

#----------------------------------------------#
# to_buyer two variables for only 12 main firms
#----------------------------------------------#
# check behavior
# temp_subsidy_type = "to_buyer"
# file_name_variable = "x1"
# size_of_subsample_temp = 0
# myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
# num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
# num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
# accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
# final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
# myests_CI_scale_X_only = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
#
# model_all_length = length(file_name_variable_list)
# myests_point_all = zeros(model_all_length,
#                          size(myests_point_scale_X_only)[1],
# 						 size(myests_point_scale_X_only)[2])
# num_correct_ineq_all = zeros(model_all_length,
#                          size(num_correct_ineq_scale_X_only)[1],
# 						 size(num_correct_ineq_scale_X_only)[2])
# num_total_ineq_all = zeros(model_all_length,
#                          size(num_total_ineq_scale_X_only)[1],
# 						 size(num_total_ineq_scale_X_only)[2])
# accuracy_all = zeros(model_all_length,
#                          size(accuracy_scale_X_only)[1])
# final_ests_point_all = zeros(model_all_length,
#                          size(final_ests_point_scale_X_only)[1])
# for ii = 1:length(file_name_variable_list)
# 	temp_subsidy_type = "to_buyer"
# 	@show temp_file_name = file_name_variable_list[ii]
# 	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
# 	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
# 	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
# 	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
# 	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(num_correct_ineq_all[ii,:,1])[2], :],
# 	                                    digits=1)
# end
#
# CI_all_table = zeros(model_all_length,2,2)
# for ii = 1:length(file_name_variable_list)
# 	temp_myests_point_all = myests_point_all[ii,num_correct_ineq_all[ii,:,1].==Int(findmax(num_correct_ineq_all[ii,:,1])[1]),:]
# 	CI_all_table[ii,:,:] = round.(hcat(
# 		Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
# 		Statistics.quantile(-temp_myests_point_all[:,2], [0.0,1.0])
# 		),digits=1)
# end
#
# total_ineq_all = Int64(num_total_ineq_all[8,findmax(accuracy_all[8,:])[2]])
# LaTeXTabulars.latex_tabular("julia_merger_table/score_results_two_variables_$(temp_subsidy_type)_subsidy_main_firms_only.tex",
# 			  Tabular("@{\\extracolsep{5pt}}lccccccccc"),
# 			  [Rule(:top),
# 			   ["","","", "", "", "", "", "", "", ""],
# 			   #["","",
# 			   #"Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
# 			   ["","",
# 			   "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]",
# 			    "[LB, UB]", "[LB, UB]", "[LB, UB]", "[LB, UB]"],
# 			   Rule(:mid),
# 			   ["Scale variables", "", "", "", "", "", "", ""],
# 			   #beta_0
# 			   ["", "", "", "", "", "", "", "", ""],
# 			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0",
# 			   "+1", "+1", "+1", "+1", "+1", "+1", "+1", "+1"],
# 			   ["" , "" , "(S)", "(S)", "(S)", "(S)",
# 			    "(S)", "(S)", "(S)", "(S)"],
# 			   #beta_1
# 			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_1",
# 			   #final_ests_point_all[1,1], "", "", "",
# 			   # "", "", "", ""],
# 			   #["" , "" ,
# 			   "[$(CI_all_table[1,1,1]), $(CI_all_table[1,2,1])]",
# 			   "", "", "",
# 			   "", "", "", ""],
# 			   #beta_2
# 			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_2",
# 				#"", final_ests_point_all[2,1], "", "",
# 				#"", "", "", ""],
# 			   #["" , "" ,
# 			   "", "[$(CI_all_table[2,1,1]), $(CI_all_table[2,2,1])]", "", "",
# 			   "", "", "", ""],
# 			   #beta_3
# 			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_3",
# 			   #"", "", final_ests_point_all[3,1], "",
# 			   #"", "", "", ""],
# 			   #["" , "" ,
# 			   "", "", "[$(CI_all_table[3,1,1]), $(CI_all_table[3,2,1])]", "",
# 			   "", "", "", ""],
# 			   # beta_4
# 			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_4",
# 			   #"", "", "", final_ests_point_all[4,1],
# 			   #"", "", "", ""],
# 			   #["" , "" ,
# 			   "", "", "", "[$(CI_all_table[4,1,1]), $(CI_all_table[4,2,1])]",
# 			   "", "", "", ""],
# 			   #beta_5
# 			   ["Share variables", "", "", "", "", "", "", "", "", ""],
# 			   ["", "", "", "", "", "", "", "", "", ""],
# 			   ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
# 			   #"", "", "", "",
# 			   #final_ests_point_all[5,1], "", "", ""],
# 			   #["" , "" ,
# 			   "", "", "", "",
# 			   "[$(CI_all_table[5,1,1]), $(CI_all_table[5,2,1])]", "", "", ""],
# 			   #beta_6
# 			   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
# 			   #"", "", "", "",
# 			   #"", final_ests_point_all[6,1], "", ""],
# 			   #["" , "" ,
# 			   "", "", "", "",
# 			   "", "[$(CI_all_table[6,1,1]), $(CI_all_table[6,2,1])]", "", ""],
# 			   #beta_7
# 			   ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
# 			   #"", "", "", "",
# 			   #"", "", final_ests_point_all[7,1], ""],
# 			   #["" , "" ,
# 			   "", "", "", "",
# 			   "", "", "[$(CI_all_table[7,1,1]), $(CI_all_table[7,2,1])]", ""],
# 			   #beta_8
# 			   ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
# 			   #"", "", "", "",
# 			   #"", "", "", final_ests_point_all[8,1]],
# 			   #["" , "" ,
# 			   "", "", "", "",
# 			   "", "", "", "[$(CI_all_table[8,1,1]), $(CI_all_table[8,2,1])]"],
# 			   #gamma merger cost (note the sign is reverse)
# 			   ["", "", "", "", "", "", "", "", "", ""],
# 			   ["", "", "", "", "", "", "", "", "", ""],
# 			   ["merger cost", "-\$\\gamma\$",
# 			   #-final_ests_point_all[1,2], -final_ests_point_all[2,2],
# 			   #-final_ests_point_all[3,2], -final_ests_point_all[4,2],
# 			   #-final_ests_point_all[5,2], -final_ests_point_all[6,2],
# 			   #-final_ests_point_all[7,2], -final_ests_point_all[8,2]],
# 			   #["" , "" ,
# 			   "[$(CI_all_table[1,1,2]), $(CI_all_table[1,2,2])]",
# 			   "[$(CI_all_table[2,1,2]), $(CI_all_table[2,2,2])]",
# 			   "[$(CI_all_table[3,1,2]), $(CI_all_table[3,2,2])]",
# 			   "[$(CI_all_table[4,1,2]), $(CI_all_table[4,2,2])]",
# 			   "[$(CI_all_table[5,1,2]), $(CI_all_table[5,2,2])]",
# 			   "[$(CI_all_table[6,1,2]), $(CI_all_table[6,2,2])]",
# 			   "[$(CI_all_table[7,1,2]), $(CI_all_table[7,2,2])]",
# 			   "[$(CI_all_table[8,1,2]), $(CI_all_table[8,2,2])]"],
# 			   #delta subsidy sensitivity
# 			   ["subsidy sensitivity", L"\delta",
# 			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
# 			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
# 			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
# 			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1],
# 			   temp_calibrated_delta_list[1], temp_calibrated_delta_list[1]],
# 			   ["", "", "", "", "", "", "", "", "", ""],
# 			   Rule(),           # a nice \hline to make it ugly
# 			   ["\$\\sharp\$ Inequalities (Point)" , "" ,
# 				total_ineq_all,
# 				total_ineq_all,
# 				total_ineq_all,
# 				total_ineq_all,
# 				total_ineq_all,
# 				total_ineq_all,
# 				total_ineq_all,
# 				total_ineq_all],
# 				#=["\\Correct Ineq" , "" ,
#  			   Int(num_correct_ineq_all[1,findmax(accuracy_all[1,:])[2]]),
#  			   Int(num_correct_ineq_all[2,findmax(accuracy_all[2,:])[2]]),
#  			   Int(num_correct_ineq_all[3,findmax(accuracy_all[3,:])[2]]),
#  			   Int(num_correct_ineq_all[4,findmax(accuracy_all[4,:])[2]]),
#  			   Int(num_correct_ineq_all[5,findmax(accuracy_all[5,:])[2]]),
#  			   Int(num_correct_ineq_all[6,findmax(accuracy_all[6,:])[2]]),
#  			   Int(num_correct_ineq_all[7,findmax(accuracy_all[7,:])[2]]),
#  			   Int(num_correct_ineq_all[8,findmax(accuracy_all[8,:])[2]])],
# 			   =#
# 			   ["\\% Inequalities" , "" ,
# 			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
# 			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
# 			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
# 			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
# 			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
# 			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
# 			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
# 			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
# 			   Rule(:bottom)])
