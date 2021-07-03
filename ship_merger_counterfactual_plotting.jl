#---------------------------------------------#
# plot for the number of groups and unmatched
#---------------------------------------------#
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

#
# @show num_group_list[1,:,:].-num_group_list[2,:,:]
# @show num_unmatched_list[1,:,:].-num_unmatched_list[2,:,:]
# @show num_unmatched_list[1,:,:].-num_unmatched_list[3,:,:]
# @show num_group_list[1,:,:].-num_group_list[3,:,:]

if length(model_list) == 3
	## parameters are set-identified so consider three scenarios
	# show num of groups
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
	# show num of unmatched
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
