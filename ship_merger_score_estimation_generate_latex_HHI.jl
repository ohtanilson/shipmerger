#---------------------#
# latex table shared
#---------------------#
# check behavior
#temp_subsidy_type = "shared"
myests_point_HHI_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_HHI_only.txt",',',Float64)
num_correct_ineq_HHI_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_HHI_only.txt",',',Float64)
num_total_ineq_HHI_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_HHI_only.txt",',',Float64)
accuracy_HHI_only = vec(num_correct_ineq_HHI_only./num_total_ineq_HHI_only)
final_ests_point_HHI_only = round.(myests_point_HHI_only[findmax(accuracy_HHI_only)[2],:],digits=2)
myests_CI_HHI_only = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_HHI_only.txt",',',Float64)
# CI
#beta = theta[1:8]
#gamma = theta[9] # coefficient on merger cost
#delta = theta[10] # coefficient on subsidy indicator
CI_all_table = zeros(2,2)
CI_all_table[:,:] = round.(hcat(
	Statistics.quantile(myests_CI_HHI_only[:,1], [0.025,0.975]),
	Statistics.quantile(-myests_CI_HHI_only[:,2], [0.025,0.975])
	),digits=1)
total_ineq_all_full_samples = 17864
total_ineq_all = Int64(num_total_ineq_HHI_only[findmax(accuracy_HHI_only[:])[2]])
LaTeXTabulars.latex_tabular("julia_merger_table/score_results_HHI_only_$(temp_subsidy_type)_subsidy.tex",
			  Tabular("@{\\extracolsep{5pt}}lcc"),
			  [Rule(:top),
			   ["","",""],
			   ["","","Point Est"],
			   ["","","[95\\% CI]"],
			   Rule(:mid),
			   #beta_0
			   ["", "", ""],
			   ["total\$_{b}\$ \$\\times\$ total\$_{t}\$", L"\beta_0","+1"],
			   ["" , "" , "(S)"],
			   #beta_1
			   ["HHI\$_{b}\$ \$\\times\$ HHI\$_{t}\$", L"\beta_{HHI}",
				final_ests_point_HHI_only[1]],
			   ["" , "" ,
			   "[$(CI_all_table[1,1]), $(CI_all_table[2,1])]",
			   "", "", "",
			   "", "", "", ""],
			   #gamma
			   ["merger cost", "-\$\\gamma\$",
			   -final_ests_point_HHI_only[2]],
			   ["" , "" ,
			   "[$(CI_all_table[1,2]), $(CI_all_table[2,2])]"],
			   #delta subsidy sensitivity
			   ["subsidy sensitivity", L"\delta",
			   #final_ests_point_scale_X_only[6], final_ests_point_scope_X_only[6], final_ests_point_full_X_only[10]],
			   temp_calibrated_delta],
			   ["", "", "", "", "", "", "", "", "", ""],
			   Rule(),           # a nice \hline to make it ugly
			   ["\$\\sharp\$ Inequalities" , "" ,
				total_ineq_all_full_samples],
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
			   round(accuracy_HHI_only[findmax(accuracy_HHI_only[:])[2]],digits=4)],
			   Rule(:bottom)])
