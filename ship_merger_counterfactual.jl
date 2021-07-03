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
using JuMP, Gurobi
using Distributions
using JLD
using LaTeXTabulars
using LaTeXStrings
#------------------------------#
#  Counterfactual
#------------------------------#
data = CSV.read("output/data_for_maximum_rank_estimation.csv",DataFrame)
rename!(data, :Column1 => :firm_id)
data_for_counterfactual = @linq data |>
	where(:type .== "(1) main")
#----------------------------------------------#
# assign output options
# It takes at least a few days for computations
# if you switch "not_run" (default) to "run".
# "not_run" justs read current restored results
# and generates Latex output for the paper.
#----------------------------------------------#
you_want_to_run_counterfactual = "not_run"
you_want_to_run_plotting = "not_run"
you_want_to_run_finding_merger_compositions = "not_run"
#-----------------------------#
# read estimated parameters
#-----------------------------#
# assign theta_hat
model_specification = "two_variables_main_firms_only" # it allows other specifications restored in julia_merger_result folder to be used
# choose subsidy specification
temp_subsidy_type = "shared"
# choose model column
file_num = 1
file_name_variable_list = ["x1","x2","x3","x4","x5","x6","x7","x8"]
file_name_variable = file_name_variable_list[file_num]
size_of_subsample_temp = 0 # only 12 main firms
calibrated_delta = 1000
myests_point_scope_X_only = readdlm("julia_merger_result/myests_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_correct_ineq_scope_X_only = readdlm("julia_merger_result/num_correct_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
num_total_ineq_scope_X_only = readdlm("julia_merger_result/num_total_ineq_subsample_size_0_$(temp_subsidy_type)_subsidy_only_$(file_name_variable)_merger_cost.txt",',',Float64)
accuracy_scope_X_only = vec(num_correct_ineq_scope_X_only./num_total_ineq_scope_X_only)
final_ests_point_scope_X_only = round.(myests_point_scope_X_only[findmax(accuracy_scope_X_only)[2],:],digits=2)
temp_myests_point_all = myests_point_scope_X_only[num_correct_ineq_scope_X_only[:,1].==Int(findmax(num_correct_ineq_scope_X_only[:,1])[1]),:]
CI_all_table = round.(hcat(
	Statistics.quantile(temp_myests_point_all[:,1], [0.0,1.0]),
	Statistics.quantile(temp_myests_point_all[:,2], [0.0,1.0])
	),digits=1)
theta_hat_all_models = round.(hcat(vcat(zeros(file_num-1),
				Statistics.quantile(temp_myests_point_all[:,1], [0.5]), # x2
				zeros(7-(file_num-1)),
				Statistics.quantile(temp_myests_point_all[:,2], [0.5]),
				calibrated_delta),
				# the cheapest scenario
				vcat(zeros(file_num-1),
				#Statistics.quantile.(myests_CI_scope_X_only[:,1:4], [0.975])[1]
				CI_all_table[2,1], # x3
				zeros(7-(file_num-1)),
				CI_all_table[1,2],
				calibrated_delta),
				# the most expensive scenario
				vcat(zeros(file_num-1),
				CI_all_table[1,1], # x3
				zeros(7-(file_num-1)),
				CI_all_table[2,2],
				calibrated_delta)
				),digits = 2)
#-----------------------------#
# set up struct and functions
#-----------------------------#
@time m = carrier_mod(randomseed = 8,
					 N = 12,
					 ϵ_dist=Distributions.Normal(0.0,5.0),
					 β_0 = 3.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 5.0,# coefficient of additional merger cost
					 ton_dim= 2)
include("ship_merger_counterfactual_functions.jl")

#---------------------------------------------------------#
# test a single computation and find a positive constant
#---------------------------------------------------------#
temp_model_id = 1
theta_hat = theta_hat_all_models[:,temp_model_id]
positive_constant = 135 # reasonable starting lower bound for all random seeds
number_of_groups_temp = 6
while number_of_groups_temp > 0
	# find positive constant term for correcting unmatched and matched
	@time utility = gen_utility_matrix_counterfactual(theta_hat,
	                        data,
							randomseed = 1,
							threshold_tonnage = 1, # 1 million
							subsidy_amount = 0,#1,#0,
							subsidy_type = "shared")
	modified_utility = copy(utility)
	# shift positive constants for capturing the unmatched environment
	for i = 1:m.N
		for j = 2:m.PN
			if IS_Sell(m,i,j)!=1
		        modified_utility[i,j] = utility[i,j] - positive_constant
			end
		end
	end
	matches = solve_equilibrium_counterfactual(m, modified_utility)
	for ii = 1:m.N
		@show ii
		@show findmax(matches[ii,:].>0)
		@show modified_utility[ii,matches[ii,:]]
	    @show m.Bundle[matches[ii,:].>0]
	end
    @show number_of_groups_temp, number_of_unmatched_temp = gen_merger_composition(m, matches)
	@show positive_constant += 1
end


#--------------------------------------------------------#
# simulate counterfactual matches for each scenario
#--------------------------------------------------------#
threshold_tonnage_list = [1, 2, 3, 4, 5, 7.5]# 6dim
subsidy_amount_list = [0, 0.1, 0.25, 0.5, 0.75, 1, 2, 3, 5] # 9dim
model_list = [1, 2, 3]
iter_end = 20
Threads.nthreads()
Threads.threadid()
JULIA_NUM_THREADS=8
println("correction positive_constant is ", positive_constant)
model_id_iter = 1
theta_hat = theta_hat_all_models[:,model_id_iter]
if you_want_to_run_counterfactual == "run"
	@time Threads.@threads for nn = 1:length(threshold_tonnage_list)
		   for mm = 1:length(subsidy_amount_list)
			number_of_groups = zeros(iter_end)
			number_of_unmatched = zeros(iter_end)
			matches_list = zeros(m.N, m.PN)
			threshold_tonnage_iter = threshold_tonnage_list[nn]
			subsidy_amount_iter = subsidy_amount_list[mm]
		    @time for kk = 1:iter_end
			    utility = gen_utility_matrix_counterfactual(theta_hat, data,
										randomseed = kk,
										threshold_tonnage = threshold_tonnage_iter, # 1 million
										subsidy_amount = subsidy_amount_iter,
										subsidy_type = "shared")
				modified_utility = copy(utility)
				#positive_constant = 135
				for i = 1:m.N
					for j = 2:m.PN
						if IS_Sell(m,i,j)!=1
					        modified_utility[i,j] = utility[i,j] - positive_constant
						end
					end
				end
			    matches = solve_equilibrium_counterfactual(m, modified_utility)
				for ii = 1:m.N
					@show ii
					@show findmax(matches[ii,:].>0)
					@show modified_utility[ii,matches[ii,:]]
				    @show m.Bundle[matches[ii,:].>0]
				end
				@show theta_hat
				@show kk, threshold_tonnage_iter, subsidy_amount_iter
			    @show number_of_groups[kk], number_of_unmatched[kk] = gen_merger_composition(m, matches)
				open("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
					DelimitedFiles.writedlm(io, convert(Array{Float64,2},matches),",")
				end
				println("\nthreshold_tonnage = $(threshold_tonnage_iter) and subsidy_amount = $(subsidy_amount_iter) iter = $(kk)\n")
				println("num of groups = $(number_of_groups[kk]) and num of unmatched = $(number_of_unmatched[kk])\n")
		    end
			open("julia_merger_result/counterfactual_number_of_groups_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
				DelimitedFiles.writedlm(io, number_of_groups,",")
			end
			open("julia_merger_result/counterfactual_number_of_unmatched_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt", "w") do io
				DelimitedFiles.writedlm(io, number_of_unmatched,",")
			end
		end
		# 6dim*9dim*1scenario*30iter on 210612
		#360402.808865 seconds (1.34 T allocations: 394.321 TiB, 33.91% gc time, 0.00% compilation time)
	end
end

# initialize result objects
kk = 1
model_id_iter = 1
threshold_tonnage_iter = 1.0
subsidy_amount_iter = 0.0
temp_d = readdlm("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
temp_iter_kk = 1
num_group_list = zeros(length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list))
num_unmatched_list = zeros(length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list))
matches_list = zeros(length(model_list),
                      length(threshold_tonnage_list),
                      length(subsidy_amount_list),
					  iter_end,
					  m.N,
					  m.PN)
#-------------------------------------#
# read simulated matching compositions
#-------------------------------------#
model_list = [1]#[1, 2, 3]
@time for model_id in model_list
	# read aggregate equilibrium numbers
	for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
		threshold_tonnage_iter = threshold_tonnage_list[nn]
		subsidy_amount_iter = subsidy_amount_list[mm]
		counterfactual_number_of_groups = readdlm("julia_merger_result/counterfactual_number_of_groups_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
		counterfactual_number_of_unmatched = readdlm("julia_merger_result/counterfactual_number_of_unmatched_model_$(model_id)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
		# take median of simulated results
		num_group_list[model_id, nn, mm] = median(counterfactual_number_of_groups)
		num_unmatched_list[model_id, nn, mm] = median(counterfactual_number_of_unmatched)
		println("model id = $(model_id) and threshold_tonnage = $(threshold_tonnage_iter) and subsidy_amount = $(subsidy_amount_iter)")
		println("num_group=",num_group_list[model_id, nn, mm],
		        ", unmatched=",num_unmatched_list[model_id, nn, mm])
	end
	# read matching compositions
	for nn = 1:length(threshold_tonnage_list), mm = 1:length(subsidy_amount_list)
		threshold_tonnage_iter = threshold_tonnage_list[nn]
		subsidy_amount_iter = subsidy_amount_list[mm]
		counterfactual_matches = zeros(m.N, m.PN, iter_end)
		for kk = 1:iter_end
	        matches_list[model_id, nn, mm, kk, :, :] = readdlm("julia_merger_result/counterfactual_matches_iter_$(kk)_model_$(model_id_iter)_threshold_tonnage_$(threshold_tonnage_iter)_subsidy_amount_$(subsidy_amount_iter).txt",',',Float64)
		end
	end
end

#---------------------------------------------#
# plot for the number of groups and unmatched
#---------------------------------------------#
if you_want_to_run_plotting == "run"
    @time include("ship_merger_counterfactual_plotting.jl")
end
#---------------------------------------------------#
# For sankey diagram
#---------------------------------------------------#
if you_want_to_run_finding_merger_compositions == "run"
    @time include("ship_merger_counterfactual_for_sankey_diagrams.jl")
end

# end script
