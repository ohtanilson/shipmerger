# this script must be located in ship_merger_score_estimation.jl
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

#want_to_run_plotting = "not_run"
@time if want_to_run_plotting == "run"
	temp_subsidy_type = "shared"
	variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
	temp_calibrated_delta_list = [50]
	# Scale variables
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
	                           calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₁", file_name_variable = "x1")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
	                           calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₂", file_name_variable = "x2")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
	                           calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₃", file_name_variable = "x3")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₄", file_name_variable = "x4")
	# Share variables
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.5:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₅", file_name_variable = "x5")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.5:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₆", file_name_variable = "x6")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.5:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₇", file_name_variable = "x7")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
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
	# Scale variables
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₁", file_name_variable = "x1")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₂", file_name_variable = "x2")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₃", file_name_variable = "x3")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-2:10:300;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₄", file_name_variable = "x4")
	# Share variables
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₅", file_name_variable = "x5")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₆", file_name_variable = "x6")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
	                           variable = "β₇", file_name_variable = "x7")
	plot_score_single_variable(temp_subsidy_type,
	                           data = data,
	                           variable_list = variable_list,
	                           domain = [-5:0.2:2;],
							   calibrated_delta_list = temp_calibrated_delta_list,
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
						  title = "Objective value: (subsidy type)=($(temp_subsidy_type)), δ=$(calibrated_delta_list)")
		end
		#Plots.vline!([0], label="", linestyle = :dash)
	end
	savefig("julia_merger_figure/plot_contour_score_two_variables_merger_cost_$(file_name_variable)_$(temp_subsidy_type)_subsidy")
end


#want_to_run_plotting = "not_run"
@time if want_to_run_plotting == "run"
	#--------------------------------#
	# shared
	#--------------------------------#
	temp_calibrated_delta_list = [20]
	#gamma_list = [70:10.0:200;]
	#temp_domain = [100.0:50.0:1500;]
	temp_subsidy_type = "shared"
	variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
	gamma_list = [1:1.0:10;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-10.0:10.0:100;],
							   #domain = [-40.0:10:40;],
							   domain = [1:20:200;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₂", file_name_variable = "x2")
	#gamma_list = [60:10.0:500;]
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
	#gamma_list = [80:10.0:500;]
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
	#gamma_list = [70:2.0:120;]
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
	#gamma_list = [380:5.0:430;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:5:15;],
							   domain = [90:2:100;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₆", file_name_variable = "x6")
	#gamma_list = [300:50.0:2000;]
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
	#gamma_list = [300:50.0:2000;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-20:10:40;],
							   domain = [-40:20:200;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₃", file_name_variable = "x3")
	#--------------------------------#
	# to_buyer
	#--------------------------------#
	temp_subsidy_type = "to_buyer"
	variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
	gamma_list = [0:1.0:10;]
	#gamma_list = [60:1.0:70;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-10.0:10.0:100;],
							   #domain = [-40.0:10:40;],
							   domain = [1:50:1000;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₂", file_name_variable = "x2")
	#gamma_list = [390:1.0:400;]
	#gamma_list = [15:5.0:200;]
	#gamma_list = [19.52:0.02:19.62;]
	#gamma_list = [400:200.0:2000;]
	#temp_domain = [1100.0:50.0:1300;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-40:10:40;],
							   domain = [1:100:1000;],
							   #domain = [-9.91:0.01:-9.8;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₄", file_name_variable = "x4")
	#gamma_list = [40.0:1.0:60.0;]
	#gamma_list = [100:50.0:300;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-20:3:1;],
							   domain = [-1000:100:0;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₅", file_name_variable = "x5")
	#gamma_list = [30.0:10.0:200.0;]
	#gamma_list = [100:10.0:130;]
	#gamma_list = [90:100.0:2000;]
	plot_contour_score_two_variables(temp_subsidy_type,
							   #domain = [-10:2:0;],
							   domain = [-1000:100:0;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₇", file_name_variable = "x7")
	#gamma_list = [40.0:10.0:200.0;]
	#gamma_list = [80:10.0:200;]
	#gamma_list = [100:50.0:300;]
	#gamma_list = [60:5.0:150;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:5:15;],
							   domain = [-1000:100:0;],
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
	#gamma_list = [36:1.0:45;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:5:15;],
							   domain = [-1000:100:0;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₆", file_name_variable = "x6")
	#gamma_list = [2.635:0.003:2.65;]
	#gamma_list = [0.3:1.0:10.0;]
	#gamma_list = [30.0:2.0:200.0;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-10:10:40;],
							   domain = [1:100:1000;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₁", file_name_variable = "x1")
	#gamma_list = [60.0:1:250;]
	#gamma_list = [300:50.0:2000;]
	plot_contour_score_two_variables(temp_subsidy_type,
	                           #domain = [-20:10:40;],
							   domain = [1:100:1000;],
							   data = data,
	                           variable_list = variable_list,
							   gamma_list = gamma_list,
							   calibrated_delta_list = temp_calibrated_delta_list,
							   variable = "β₃", file_name_variable = "x3")
end
#3039.685852 seconds (11.68 G allocations: 1.109 TiB, 4.31% gc time)
