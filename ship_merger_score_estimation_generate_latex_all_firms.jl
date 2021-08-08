#---------------------#
# latex table shared
#---------------------#
# check behavior
#temp_subsidy_type = "shared"
file_name_variable = "x1"
size_of_subsample_temp = 106#30
myests_point_scale_X_only = readdlm("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_correct_ineq_scale_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_total_ineq_scale_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
accuracy_scale_X_only = vec(num_correct_ineq_scale_X_only./num_total_ineq_scale_X_only)
final_ests_point_scale_X_only = round.(myests_point_scale_X_only[findmax(accuracy_scale_X_only)[2],:],digits=2)
myests_CI_scale_X_only = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)

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
	#temp_subsidy_type = "shared"
	@show temp_file_name = file_name_variable_list[ii]
	myests_point_all[ii,:,:] = readdlm("julia_merger_result/myests_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_correct_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_correct_ineq_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	num_total_ineq_all[ii,:,:] = readdlm("julia_merger_result/num_total_ineq_subsample_size_$(size_of_fullsample)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
	accuracy_all[ii,:] = num_correct_ineq_all[ii,:,1]./num_total_ineq_all[ii,:,1]
	final_ests_point_all[ii,:] = round.(myests_point_all[ii, findmax(accuracy_all[ii,:])[2], :],
	                                    digits=2)
	myests_CI_all[ii,:,:] = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_only_$(temp_file_name)_merger_cost.txt",',',Float64)
end
maximum(num_correct_ineq_all[:,:,1],dims = 2)
final_ests_point_all
total_ineq_all_full_samples = 17864
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
               ["","","(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)"],
			   ["","",
			   "Point Est", "Point", "Point", "Point", "Point", "Point", "Point", "Point"],
			   ["","",
			   "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]",
			    "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]"],
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
			   ["Share variables", "", "", "", "", "", "", "", "", ""],
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
				total_ineq_all_full_samples,
				total_ineq_all_full_samples,
				total_ineq_all_full_samples,
				total_ineq_all_full_samples,
				total_ineq_all_full_samples,
				total_ineq_all_full_samples,
				total_ineq_all_full_samples,
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
			   round(accuracy_all[1,findmax(accuracy_all[1,:])[2]],digits=4),
			   round(accuracy_all[2,findmax(accuracy_all[2,:])[2]],digits=4),
			   round(accuracy_all[3,findmax(accuracy_all[3,:])[2]],digits=4),
			   round(accuracy_all[4,findmax(accuracy_all[4,:])[2]],digits=4),
			   round(accuracy_all[5,findmax(accuracy_all[5,:])[2]],digits=4),
			   round(accuracy_all[6,findmax(accuracy_all[6,:])[2]],digits=4),
			   round(accuracy_all[7,findmax(accuracy_all[7,:])[2]],digits=4),
			   round(accuracy_all[8,findmax(accuracy_all[8,:])[2]],digits=4)],
			   Rule(:bottom)])

sample_production_level = round.(
              [mean_of_x0_from_summary_stats^2 + final_ests_point_all[1,1]*mean_of_x_from_summary_stats[1]^2 - final_ests_point_all[1,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[2,1]*mean_of_x_from_summary_stats[2]^2 - final_ests_point_all[2,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[3,1]*mean_of_x_from_summary_stats[3]^2 - final_ests_point_all[3,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[4,1]*mean_of_x_from_summary_stats[4]^2 - final_ests_point_all[4,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[5,1]*mean_of_x_from_summary_stats[5]^2 - final_ests_point_all[5,2],
 			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[6,1]*mean_of_x_from_summary_stats[6]^2 - final_ests_point_all[6,2],
 			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[7,1]*mean_of_x_from_summary_stats[7]^2 - final_ests_point_all[7,2],
 			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[8,1]*mean_of_x_from_summary_stats[8]^2 - final_ests_point_all[8,2]],
			   digits = 3)

sample_production_level_UP_1SD = round.(
              [mean_of_x0_from_summary_stats^2 + final_ests_point_all[1,1]*(max(mean_of_x_from_summary_stats[1] + sd_of_x_from_summary_stats[1],0))^2 - final_ests_point_all[1,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[2,1]*(max(mean_of_x_from_summary_stats[2] + sd_of_x_from_summary_stats[2],0))^2 - final_ests_point_all[2,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[3,1]*(max(mean_of_x_from_summary_stats[3] + sd_of_x_from_summary_stats[3],0))^2 - final_ests_point_all[3,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[4,1]*(max(mean_of_x_from_summary_stats[4] + sd_of_x_from_summary_stats[4],0))^2 - final_ests_point_all[4,2],
			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[5,1]*(max(mean_of_x_from_summary_stats[5] + sd_of_x_from_summary_stats[5],0))^2 - final_ests_point_all[5,2],
 			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[6,1]*(max(mean_of_x_from_summary_stats[6] + sd_of_x_from_summary_stats[6],0))^2 - final_ests_point_all[6,2],
 			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[7,1]*(max(mean_of_x_from_summary_stats[7] + sd_of_x_from_summary_stats[7],0))^2 - final_ests_point_all[7,2],
 			   mean_of_x0_from_summary_stats^2 + final_ests_point_all[8,1]*(max(mean_of_x_from_summary_stats[8] + sd_of_x_from_summary_stats[8],0))^2 - final_ests_point_all[8,2]],
			   digits = 3)
sample_production_level_diff = round.(sample_production_level_UP_1SD.-sample_production_level,digits = 3)
LaTeXTabulars.latex_tabular("julia_merger_table/ratio_score_results_two_variables_$(temp_subsidy_type)_subsidy.tex",
			  Tabular("@{\\extracolsep{5pt}}lcccccccc"),
			  [Rule(:top),
			   ["\$m\$",
			   "\$x_0\$",
			   "\$\\beta_m\$",
			   "\$\\bar{x}_m\$",
			   "\$1SD(x_m)\$",
			   "\$-\\gamma\$",
			   "\$V(\\bar{x}_m)=x_0 x_0 + \\beta_m \\bar{x}_m \\bar{x}_m-\\gamma\$",
			   "\$V(\\bar{x}_m+1SD(x_m))\$",
			   "\$V(\\bar{x}_m+1SD(x_m))-V(\\bar{x}_m)\$"],
			   Rule(:mid),
			   ["1", mean_of_x0_from_summary_stats, final_ests_point_all[1,1], mean_of_x_from_summary_stats[1], sd_of_x_from_summary_stats[1],
			    -final_ests_point_all[1,2], sample_production_level[1], sample_production_level_UP_1SD[1], sample_production_level_diff[1]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["2", mean_of_x0_from_summary_stats, final_ests_point_all[2,1], mean_of_x_from_summary_stats[2], sd_of_x_from_summary_stats[2],
			    -final_ests_point_all[2,2], sample_production_level[2], sample_production_level_UP_1SD[2], sample_production_level_diff[2]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["3", mean_of_x0_from_summary_stats, final_ests_point_all[3,1], mean_of_x_from_summary_stats[3], sd_of_x_from_summary_stats[3],
			    -final_ests_point_all[3,2], sample_production_level[3], sample_production_level_UP_1SD[3], sample_production_level_diff[3]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["4", mean_of_x0_from_summary_stats, final_ests_point_all[4,1], mean_of_x_from_summary_stats[4], sd_of_x_from_summary_stats[4],
			    -final_ests_point_all[4,2], sample_production_level[4], sample_production_level_UP_1SD[4], sample_production_level_diff[4]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["", "", "", "", "", "", "", "", ""],
			   ["5", mean_of_x0_from_summary_stats, final_ests_point_all[5,1], mean_of_x_from_summary_stats[5], sd_of_x_from_summary_stats[5],
			   -final_ests_point_all[5,2], sample_production_level[5], sample_production_level_UP_1SD[5], sample_production_level_diff[5]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["6", mean_of_x0_from_summary_stats, final_ests_point_all[6,1], mean_of_x_from_summary_stats[6], sd_of_x_from_summary_stats[6],
			   -final_ests_point_all[6,2], sample_production_level[6], sample_production_level_UP_1SD[6], sample_production_level_diff[6]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["7", mean_of_x0_from_summary_stats, final_ests_point_all[7,1], mean_of_x_from_summary_stats[7], sd_of_x_from_summary_stats[7],
			   -final_ests_point_all[7,2], sample_production_level[7], sample_production_level_UP_1SD[7], sample_production_level_diff[7]],
			   ["", "", "", "", "", "", "", "", ""],
			   ["8", mean_of_x0_from_summary_stats, final_ests_point_all[8,1], mean_of_x_from_summary_stats[8], sd_of_x_from_summary_stats[8],
			    -final_ests_point_all[8,2], sample_production_level[8], sample_production_level_UP_1SD[8], sample_production_level_diff[8]],
			   ["", "", "", "", "", "", "", "", ""],
			   Rule(:bottom)])
