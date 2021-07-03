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
# setup module
#-----------------------#
@time m = carrier_mod(randomseed = 8,
					 N = 8,
					 β_0 = 3.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 5.0,# coefficient of additional merger cost
					 ton_dim= 2)

#------------------------------#
#  load data and functions
#------------------------------#
data = CSV.read("output/data_for_maximum_rank_estimation.csv",DataFrame)
rename!(data, :Column1 => :firm_id)
gdf = groupby(data, :group)
temp_info_sum = Dict([("liner_sum", combine(gdf, [:liner] => sum)),
                 ("special_sum", combine(gdf, [:special] => sum)),
				 ("tanker_sum", combine(gdf, [:tanker] => sum)),
				 ("tramper_sum", combine(gdf, [:tramper] => sum))])
include("ship_merger_score_estimation_point_estimates_and_confidence_intervals.jl")

# assign subsidy type
temp_subsidy_type = "shared" # or "to_buyer"

#----------------------------------------------#
# assign output options
# It takes at least a few days for each estimation block
# if you switch "not_run" (default) to "run".
# "not_run" justs read current restored results
# and generates Latex output for the paper.
#----------------------------------------------#
want_to_run_plotting = "not_run" # change to "run" if you want to generate all plots
want_to_run_multivariate_estimation = "not_run" # change to "run" if you want to run estimation
want_to_run_two_variable_estimation = "not_run"
want_to_run_only_main_firms = "not_run"

#----------------------------#
# Plotting objective function
#----------------------------#
include("ship_merger_score_estimation_plotting.jl")
#------------------------------------------#
# point estimate of multivariate estimation
#------------------------------------------#
size_of_fullsample = length(data.firm_id)-12
common_num_its_bootstrap = 200
common_num_steps_DE = 100

@time if want_to_run_multivariate_estimation == "run"
	@time point_estimate(temp_subsidy_type,
	                     num_steps_DE_temp = common_num_steps_DE,#50,
	                     num_its_temp = 50,#20,
						 #num_its_temp = 2,
						 size_of_fullsample = size_of_fullsample,
						 temp_temp_calibrated_delta = 5)
	#@time point_estimate(num_steps_DE_temp = 100, num_its_temp = 2)
	#1007.291976 seconds (5.74 G allocations: 557.036 GiB, 6.62% gc time)
	#@time point_estimate(num_steps_DE_temp = 100,  num_its_temp = 10)
	#3812.250026 seconds (20.37 G allocations: 2.009 TiB, 6.69% gc time)
	#@time point_estimate(num_steps_DE_temp = 50,  num_its_temp = 100)
	#151653.998867 seconds (160.68 G allocations: 15.258 TiB, 6.47% gc time)
end
@time if want_to_run_multivariate_estimation == "run"
	@time construct_CI(temp_subsidy_type,
	                   num_its_bootstrap = common_num_its_bootstrap,
	                   #num_its_bootstrap = 2,
	                   num_steps_DE_temp = common_num_steps_DE,
					   size_of_subsample_temp = size_of_fullsample,
					   temp_temp_calibrated_delta = 5,
					   temp_sampling_method = "bootstrap")
	#34558.137605 seconds (162.34 G allocations: 15.437 TiB, 5.96% gc time)
	#@time point_estimate(num_steps_DE_temp = 50,  num_its_bootstrap = 200)
	#73597.442566 seconds (90.89 G allocations: 8.674 TiB, 11.04% gc time)
end
# generate latex output
include("ship_merger_score_estimation_generate_latex_multivariates.jl")

#-----------------------#
# estimatetwo_variables
#-----------------------#
variable_list = ["β₁","β₂","β₃","β₄","β₅","β₆","β₇","β₈","γ","δ"]
file_name_variable_list = ["x1","x2","x3","x4","x5","x6","x7","x8"]
Threads.nthreads()
JULIA_NUM_THREADS=8
temp_calibrated_delta_list = [5]
@time if want_to_run_two_variable_estimation == "run"
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
	    @time point_estimate_two_variables(temp_subsidy_type;
							data = data,
	    					variable_list = variable_list,
		                    size_of_fullsample = size_of_fullsample,
		                    num_steps_DE_temp = common_num_steps_DE,
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

@time if want_to_run_two_variable_estimation == "run"
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
		construct_CI_two_variables(temp_subsidy_type;
							  data = data,
    	    				  variable_list = variable_list,
			                  num_its_bootstrap = common_num_its_bootstrap,#100,
			                  num_steps_DE_temp = common_num_steps_DE,
			                  size_of_subsample_temp = size_of_fullsample,
							  calibrated_delta_list = temp_calibrated_delta_list,
							  variable = temp_target_variable,
							  file_name_variable = temp_file_name,
  							  info_sum = temp_info_sum,
							  temp_sampling_method = "bootstrap")
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
# generate latex output
include("ship_merger_score_estimation_generate_latex_all_firms.jl")

#-----------------------#
# estimate only main firm
#-----------------------#
temp_calibrated_delta_list = [1000]
@time if want_to_run_only_main_firms == "run"
	# choose only main firms
	@time Threads.@threads for ii = 1:length(file_name_variable_list)
		temp_file_name = file_name_variable_list[ii]
		temp_target_variable = variable_list[ii]
	    @time point_estimate_two_variables(temp_subsidy_type;
		                    data = data,
		                    variable_list = variable_list,
		                    size_of_fullsample = 0, # so we have only 12 main firms
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
# generate latex output
include("ship_merger_score_estimation_generate_latex_only_main_firms.jl")
