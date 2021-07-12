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
include("ship_merger_score_estimation_point_estimates_and_confidence_intervals_HHI.jl")

# assign subsidy type
temp_subsidy_type = "shared" # or "to_buyer"

#----------------------------------------------#
# assign output options
# It takes at least a few days for each estimation block
# if you switch "not_run" (default) to "run".
# "not_run" justs read current restored results
# and generates Latex output for the paper.
#----------------------------------------------#
want_to_run_only_HHI_estimation = "run" # change to "run" if you want to run estimation

#------------------------------------------#
# point estimate of multivariate estimation
#------------------------------------------#
size_of_fullsample = length(data.firm_id)-12
common_num_its_bootstrap = 200
common_num_steps_DE = 100
temp_calibrated_delta = 5
@time if want_to_run_only_HHI_estimation == "run"
	@time point_estimate_only_HHI(temp_subsidy_type,
	                     num_steps_DE_temp = common_num_steps_DE,#50,
	                     num_its_temp = 50,#20,
						 #num_its_temp = 2,
						 size_of_fullsample = size_of_fullsample,
						 temp_temp_calibrated_delta = temp_calibrated_delta)
	#@time point_estimate(num_steps_DE_temp = 100, num_its_temp = 50)
	#73597.442566 seconds (90.89 G allocations: 8.674 TiB, 11.04% gc time)
end
@time if want_to_run_only_HHI_estimation == "run"
	@time construct_CI_only_HHI(temp_subsidy_type,
	                   num_its_bootstrap = common_num_its_bootstrap,
	                   #num_its_bootstrap = 2,
	                   num_steps_DE_temp = common_num_steps_DE,
					   size_of_subsample_temp = size_of_fullsample,
					   temp_temp_calibrated_delta = temp_calibrated_delta,
					   temp_sampling_method = "bootstrap")
	#34558.137605 seconds (162.34 G allocations: 15.437 TiB, 5.96% gc time)
	#@time point_estimate(num_steps_DE_temp = 50,  num_its_bootstrap = 200)
	#63609.682311 seconds (276.74 G allocations: 18.689 TiB, 4.73% gc time, 0.00% compilation time)
end
# generate latex output
include("ship_merger_score_estimation_generate_latex_HHI.jl")
