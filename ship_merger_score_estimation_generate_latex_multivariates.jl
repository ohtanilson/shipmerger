# this script must be located in ship_merger_score_estimation.jl
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
	myests_CI_scale_X_only = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scale_X_only.txt",',',Float64)
	myests_CI_scope_X_only = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_scope_X_only.txt",',',Float64)
	myests_CI_full_X_only = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_full_X_only.txt",',',Float64)
	myests_CI_x45678_merger_cost = readdlm("julia_merger_result/CI_myests_subsample_size_$(size_of_subsample_temp)_$(temp_subsidy_type)_subsidy_x45678_merger_cost.txt",',',Float64)

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
	LaTeXTabulars.latex_tabular("julia_merger_table/score_results_multivariate_$(temp_subsidy_type)_subsidy.tex",
	              Tabular("@{\\extracolsep{5pt}}lccccc"),
	              [Rule(:top),
				   ["","","", "Value Function", "", ""],
	               ["","",
				   "Point Estimate", "Point Estimate", "Point Estimate", "Point Estimate"],
				   ["","",
				   "[95\\% CI]", "[95\\% CI]", "[95\\% CI]", "[95\\% CI]"],
				   Rule(:mid),

				   ["Scale variables", "", "", "", "", ""],
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
				   ["Share variables", "", "", "", "", ""],
				   ["", "", "", "", "", ""],
	               ["liner\$_{b}\$ \$\\times\$ liner\$_{t}\$", L"\beta_5",
				    "", final_ests_point_scope_X_only[1], final_ests_point_full_X_only[5],
					 final_ests_point_x45678_merger_cost[2]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,1]), $(CI_scope_X_only_table[2,1])]",
				   "[$(CI_full_X_only_table[1,5]), $(CI_full_X_only_table[2,5])]",
				   "[$(CI_x45678_merger_cost_table[1,2]), $(CI_x45678_merger_cost_table[2,2])]"],
				   #beta_6
				   ["tramper\$_{b}\$ \$\\times\$ tramper\$_{t}\$", L"\beta_6",
				    "", final_ests_point_scope_X_only[2], final_ests_point_full_X_only[6],
					 final_ests_point_x45678_merger_cost[3]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,2]), $(CI_scope_X_only_table[2,2])]",
				   "[$(CI_full_X_only_table[1,6]), $(CI_full_X_only_table[2,6])]",
				   "[$(CI_x45678_merger_cost_table[1,3]), $(CI_x45678_merger_cost_table[2,3])]"],
				   #beta_7
	               ["special\$_{b}\$ \$\\times\$ special\$_{t}\$", L"\beta_7",
				    "", final_ests_point_scope_X_only[3], final_ests_point_full_X_only[7],
					 final_ests_point_x45678_merger_cost[4]],
				   ["" , "" ,
				   "",
				   "[$(CI_scope_X_only_table[1,3]), $(CI_scope_X_only_table[2,3])]",
				   "[$(CI_full_X_only_table[1,7]), $(CI_full_X_only_table[2,7])]",
				   "[$(CI_x45678_merger_cost_table[1,4]), $(CI_x45678_merger_cost_table[2,4])]"],
				   #beta_8
	               ["tanker\$_{b}\$ \$\\times\$ tanker\$_{t}\$", L"\beta_8",
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
	               ["subsidy sensitivity", L"\delta",
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

# model 1234
generate_score_table_model_1234(temp_subsidy_type = temp_subsidy_type,
								size_of_subsample_temp = size_of_fullsample)
# find point-estimate LB of model 2
if 1 == 2 # if you want to check, switch 1 to 2
    @time plot_point_estimate_histogram(temp_subsidy_type)
	#0.201591 seconds (261.57 k allocations: 10.609 MiB)
    # @time plot_point_estimate_histogram("to_buyer")
end
