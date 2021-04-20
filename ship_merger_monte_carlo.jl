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
using Distributions
using Dates
include(joinpath(dirname(@__FILE__),"functions.jl"))
using .functions
temp_randomseed = 1
#--------------------------#
# Identification plots
#--------------------------#
@time m = carrier_mod(randomseed = 3,#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)
@show m.ton
temp_threshold_tonnage = 1
temp_subsidy_amount = 1
temp_subsidy_type = "shared"
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = temp_subsidy_type)
utility_to_buyer, obsd_to_buyer = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = "to_buyer")
@show obsd.tarid
@show obsd_to_buyer.tarid

m.Bundle[obsd_to_buyer.tarid]
#--------------------------#
# (1) s_shared
#--------------------------#
temp_threshold_tonnage = 1
temp_subsidy_amount = 1
temp_subsidy_type = "shared"
@time m = carrier_mod(randomseed = 3,#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = temp_subsidy_type)
m.Bundle[obsd.tarid]
obsd.total_size_target + obsd.total_size_buyer
total_unmatched = sum(obsd.tarid .== 1)
total_num_group = length(unique(obsd.tarid[obsd.tarid .!= 1]))
domain = [-10:1.0:40;]
res_β_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_β_0[i], ~ = score_b(m, [iter, m.δ_0, m.γ_0], obsd,
	                m.N, temp_threshold_tonnage,
					temp_subsidy_amount,
					subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_β_0,
           label="score",
           xlabel="β where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
           title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.vline!([m.β_0], label="true where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]", linestyle = :dot)
savefig("julia_merger_figure/plot_beta_$(temp_subsidy_type)_subsidy")

domain = [-10:1.0:40;]
res_δ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_δ_0[i], ~ = score_b(m, [m.β_0, iter, m.γ_0], obsd,
	                        m.N, temp_threshold_tonnage,
						    temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_δ_0, label="score",
           xlabel="δ (subsidy sensitivity) where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
           title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type)) ")
Plots.vline!([m.δ_0], label="true where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]", linestyle = :dot)
savefig("julia_merger_figure/plot_delta_$(temp_subsidy_type)_subsidy")

domain = [-10:1.0:40;]
res_γ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_γ_0[i], ~ = score_b(m, [m.β_0, m.δ_0, iter],
	                        obsd, m.N, temp_threshold_tonnage,
							temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_γ_0, label="score",
           xlabel = "γ (merger cost) where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
           title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type)) ")
Plots.vline!([m.γ_0], label="true where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]", linestyle = :dot)
savefig("julia_merger_figure/plot_gamma_$(temp_subsidy_type)_subsidy")

domain1 = [-10:1.0:30;]
domain2 = [-10:1.0:30;]
res_contor = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor[i,j], ~ = score_b(m, [domain1[i], m.δ_0, domain2[j]],
	                             obsd, m.N, temp_threshold_tonnage,
								 temp_subsidy_amount,
								 subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor', fill = true,
              ylabel = "γ (merger cost)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)], $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.hline!([m.γ_0], label="true γ", linestyle = :dash)
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_gamma_$(temp_subsidy_type)_subsidy")

domain1 = [-10:1.0:30;]
domain2 = [-10:1.0:30;]
res_contor2 = zeros(length(domain1),length(domain2))
@time for i in 1:length(domain1),j in 1:length(domain2)
	res_contor2[i,j], ~ = score_b(m, [m.β_0, domain1[i], domain2[j]],
	                              obsd, m.N, temp_threshold_tonnage,
								  temp_subsidy_amount,
								  subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor2', fill = true,
              ylabel = "γ (merger cost)",
			  xlabel = "δ (subsidy sensitivity) where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash)
savefig("julia_merger_figure/contour_gamma_delta_$(temp_subsidy_type)_subsidy")

domain1 = [-10:1.0:30;]
domain2 = [-10:1.0:30;]
res_contor3 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor3[i,j], ~ = score_b(m, [domain1[i], domain2[j], m.γ_0],
	                              obsd, m.N, temp_threshold_tonnage,
								  temp_subsidy_amount,
								  subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor3', fill = true,
              ylabel = "δ (subsidy sensitivity)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)], $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_delta_$(temp_subsidy_type)_subsidy")

#--------------------------#
# (2) s_to_buyer
#--------------------------#
temp_subsidy_type = "to_buyer"
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = temp_subsidy_type)
m.Bundle[obsd.tarid]
total_unmatched = sum(obsd.tarid .== 1)
total_num_group = length(unique(obsd.tarid[obsd.tarid .!= 1]))
domain = [-10:5.0:40;]
res_β_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_β_0[i], ~ = score_b(m, [iter, m.δ_0, m.γ_0], obsd,
	                m.N, temp_threshold_tonnage,
					temp_subsidy_amount,
					subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_β_0,
           label="score",
           xlabel="β where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
           title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.vline!([m.β_0], label="true where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]", linestyle = :dot)
savefig("julia_merger_figure/plot_beta_$(temp_subsidy_type)_subsidy")

domain = [-10:1.0:40;]
res_δ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_δ_0[i], ~ = score_b(m, [m.β_0, iter, m.γ_0], obsd,
	                        m.N, temp_threshold_tonnage,
						    temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_δ_0, label="score",
           xlabel="δ (subsidy sensitivity) where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
           title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type)) ")
Plots.vline!([m.δ_0], label="true where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]", linestyle = :dot)
savefig("julia_merger_figure/plot_delta_$(temp_subsidy_type)_subsidy")

domain = [-10:1.0:40;]
res_γ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_γ_0[i], ~ = score_b(m, [m.β_0, m.δ_0, iter],
	                        obsd, m.N, temp_threshold_tonnage,
							temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_γ_0, label="score",
           xlabel = "γ (merger cost) where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
           title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type)) ")
Plots.vline!([m.γ_0], label="true where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]", linestyle = :dot)
savefig("julia_merger_figure/plot_gamma_$(temp_subsidy_type)_subsidy")

domain1 = [-10:1.0:40;]
domain2 = [-10:1.0:40;]
res_contor = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor[i,j], ~ = score_b(m, [domain1[i], m.δ_0, domain2[j]],
	                             obsd, m.N, temp_threshold_tonnage,
								 temp_subsidy_amount,
								 subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor', fill = true,
              ylabel = "γ (merger cost)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)], $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.hline!([m.γ_0], label="true γ", linestyle = :dash)
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_gamma_$(temp_subsidy_type)_subsidy")

domain1 = [-10:1.0:40;]
domain2 = [-10:1.0:40;]
res_contor2 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor2[i,j], ~ = score_b(m, [m.β_0, domain1[i], domain2[j]],
	                              obsd, m.N, temp_threshold_tonnage,
								  temp_subsidy_amount,
								  subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor2', fill = true,
              ylabel = "γ (merger cost)",
			  xlabel = "δ (subsidy sensitivity) where $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash)
savefig("julia_merger_figure/contour_gamma_delta_$(temp_subsidy_type)_subsidy")

domain1 = [-10:1.0:40;]
domain2 = [-10:1.0:40;]
res_contor3 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor3[i,j], ~ = score_b(m, [domain1[i], domain2[j], m.γ_0],
	                              obsd, m.N, temp_threshold_tonnage,
								  temp_subsidy_amount,
								  subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor3', fill = true,
              ylabel = "δ (subsidy sensitivity)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)], $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_delta_$(temp_subsidy_type)_subsidy")

#--------------------------#
# (3) density of estimates
#--------------------------#
function maxscore_mc_temp(num_agents::Int64;
                     num_its::Int64 = 10,
                     threshold_tonnage = 100,
					 subsidy_amount = 10,
					 subsidy_type = "to_buyer",
					 β_0_temp = 0.0,
					 δ_0_temp = 1.0,
					 γ_0_temp = 4.0,
					 search_range_temp = (-20.0, 20.0))
    start = Dates.unix2datetime(time())
	param_dim = 3
    myests = zeros(num_its, param_dim)
	mynum_correct_ineq = zeros(num_its, 1)
	truenum_correct_ineq = zeros(num_its, 1)
	num_all_ineq = zeros(num_its)
	solved_case_index = zeros(Int64, 0)
	subsidy_used_or_not_index = zeros(num_its)
	all_unmatched_or_not_index = zeros(num_its)
    for i = 1:num_its
	   @show i
       println("Create obsdat for iteration $i \n" )
	   @time m = carrier_mod(randomseed = i,
							N = num_agents,
							β_0 = β_0_temp,# coefficient of covariates
	  					    δ_0 = δ_0_temp,# coefficient of subsidy
	  					    γ_0 = γ_0_temp,# coefficient of additional merger cost
							ton_dim= 2)
	   #matches = solve_equilibrium(m,gen_utility_matrix)
	   #utility, matches = solve_equilibrium(m, threshold_tonnage)
	   utility, obsd = gen_data(m, threshold_tonnage = threshold_tonnage,
	                            subsidy_amount = subsidy_amount,
							    subsidy_type = subsidy_type)
	   @show matching_index = hcat(obsd.buyid, obsd.tarid)
	   if size(obsd)[1] > m.N
		   	println("skip because the matching outcome is not integer")
	   else
		   if sum(obsd.total_size_target + obsd.total_size_buyer .>1) > 0
			   println("simulated data does not use subsidy")
			   subsidy_used_or_not_index[i] = 1
		   end
		   if sum(obsd.tarid) == m.N
			   println("simulated data observe no mergers")
			   all_unmatched_or_not_index[i] = 1
		   end
		   solved_case_index = vcat(solved_case_index, i)
		   function score_bthis(beta::Vector{Float64}, obsdat::DataFrame)
			   score, ~, ~ = score_b(m, beta, obsdat, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
			   res = -1.0*score + 100000.0 # need to be Float64 for bboptimize
			   #println("score:$(res) \n" )
			   return res
		   end
		   m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsd);
		                                    SearchRange = search_range_temp,
		                                    NumDimensions = param_dim,
											Method = :de_rand_1_bin,
											MaxSteps = 50)
		   println("score: ", score_bthis(m_res.archive_output.best_candidate, obsd))
		   mynum_correct_ineq[i,1],~ = score_b(m, m_res.archive_output.best_candidate, obsd, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
		   truenum_correct_ineq[i,1],total_num_ineq = score_b(m, [m.β_0, m.δ_0, m.γ_0], obsd, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
		   println("the number of correct inequalities: ", mynum_correct_ineq[i,1])
		   println("the TRUE number of score of correct: ", truenum_correct_ineq[i,1])
		   myests[i,:] = m_res.archive_output.best_candidate
		   num_all_ineq[i] = total_num_ineq
	   end
   end
   #runtime = proc.time() - ptm
   finish = convert(Int, Dates.value(Dates.unix2datetime(time())-start))/1000
   #mean = mean(vec-meanmat)
   #sqrt = sqrt(mean((vec-meanmat)^2))
   meanmat = reshape(repeat([m.β_0, m.δ_0, m.γ_0], num_its), param_dim, num_its)' # true parameter
   res_abs_mean_err = mean(abs.((myests.-meanmat)[solved_case_index,:]), dims = 1)
   res_sqrt = mean(((myests.-meanmat)[solved_case_index,:]).^2, dims = 1)
   mynum_correct_ineq = mynum_correct_ineq[solved_case_index,:]
   truenum_correct_ineq = truenum_correct_ineq[solved_case_index,:]
   num_all_ineq = num_all_ineq[solved_case_index]
   subsidy_used_or_not_index
   #@show myests
   return myests, res_abs_mean_err, res_sqrt, mynum_correct_ineq, truenum_correct_ineq, num_all_ineq, subsidy_used_or_not_index, all_unmatched_or_not_index
end


@time m = carrier_mod(randomseed = 1,
					 N = 8,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)

temp_simulate_num = 1000
#----------------------------------------------#
# plot subsidy (to buyer)
#----------------------------------------------#
you_want_run = "not_run"
if you_want_run == "run"
    @time res_theta, ~, ~, ~, ~, ~, subsidy_used_or_not_index, all_unmatched_or_not_index = maxscore_mc_temp(8,num_its = temp_simulate_num,
                                                                                                               threshold_tonnage = 1,
																											   subsidy_amount = 1,
																											   subsidy_type = "to_buyer")
	open("julia_merger_result/res_theta_to_buyer.txt", "w") do io
		DelimitedFiles.writedlm(io, res_theta,",")
	end
	open("julia_merger_result/subsidy_used_or_not_index_to_buyer.txt", "w") do io
		DelimitedFiles.writedlm(io, subsidy_used_or_not_index,",")
	end
end
res_theta = readdlm("julia_merger_result/res_theta_to_buyer.txt",',',Float64)
subsidy_used_or_not_index = vec(readdlm("julia_merger_result/subsidy_used_or_not_index_to_buyer.txt",',',Float64))
#res_theta_solved = res_theta[res_theta.!=0]
res_theta_solved = copy(res_theta)
for i = 1:size(res_theta)[1], j = 1:size(res_theta)[2]
	if res_theta[i,j] == 0.0
		res_theta_solved[i,j] = copy(NaN)
	end
end
Plots.histogram(res_theta_solved[subsidy_used_or_not_index.<1,1],
                bins = 60,
                ylabel = "count",
                xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated β",
				title = "subsidy (to buyer)",
				xlim = [-20,20],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[subsidy_used_or_not_index.>0,1],
                bins = 60,
                label = "estimated β (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.β_0], label="true β", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_beta_to_buyer")

Plots.histogram(res_theta_solved[subsidy_used_or_not_index.<1,2],
                bins = 60,
                ylabel = "count",
				xlabel = "δ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
				label = "estimated δ",
				title = "subsidy (to buyer)",
				xlim = [-20,20],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[subsidy_used_or_not_index.>0,2],
                bins = 60,
                label = "estimated δ (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.δ_0], label="true δ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_delta_to_buyer")

Plots.histogram(res_theta_solved[subsidy_used_or_not_index.<1,3],
                bins = 60,
                ylabel = "count",
				xlabel = "γ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated γ",
				title = "subsidy (to buyer)",
				xlim = [-20,20],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[subsidy_used_or_not_index.>0,3],
                bins = 60,
                label = "estimated γ (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_gamma_to_buyer")
# Create Table
nanmedian(x) = median(filter(!isnan,x))
nanmean(x) = mean(filter(!isnan,x))
RMSE(x,true_param) = sqrt(nanmean((x.-true_param).^2))
median_bias_beta_with_subsidy_to_buyer = nanmedian(res_theta_solved[subsidy_used_or_not_index.<1,1].-m.β_0)
median_bias_beta_without_subsidy_to_buyer = nanmedian(res_theta_solved[subsidy_used_or_not_index.>0,1].-m.β_0)
median_bias_gamma_with_subsidy_to_buyer = nanmedian(res_theta_solved[subsidy_used_or_not_index.<1,3].-m.γ_0)
median_bias_gamma_without_subsidy_to_buyer = nanmedian(res_theta_solved[subsidy_used_or_not_index.>0,3].-m.γ_0)
RMSE_beta_with_subsidy_to_buyer = RMSE(res_theta_solved[subsidy_used_or_not_index.<1,1],m.β_0)
RMSE_beta_without_subsidy_to_buyer = RMSE(res_theta_solved[subsidy_used_or_not_index.>0,1],m.β_0)
RMSE_gamma_with_subsidy_to_buyer = RMSE(res_theta_solved[subsidy_used_or_not_index.<1,3],m.γ_0)
RMSE_gamma_without_subsidy_to_buyer = RMSE(res_theta_solved[subsidy_used_or_not_index.>0,3],m.γ_0)

#----------------------------------------------#
# plot subsidy (shared)
#----------------------------------------------#
you_want_run = "not_run"
if you_want_run == "run"
	@time res_theta, ~, ~, ~, ~, ~, subsidy_used_or_not_index, all_unmatched_or_not_index = maxscore_mc_temp(8, num_its = temp_simulate_num,
	                                                                                                                      threshold_tonnage = 100,
																														  subsidy_amount = 10,
																														  subsidy_type = "shared")
	open("julia_merger_result/res_theta_shared.txt", "w") do io
		DelimitedFiles.writedlm(io, res_theta,",")
	end
	open("julia_merger_result/subsidy_used_or_not_index_shared.txt", "w") do io
		DelimitedFiles.writedlm(io, subsidy_used_or_not_index,",")
	end
end

res_theta = readdlm("julia_merger_result/res_theta_shared.txt",',',Float64)
subsidy_used_or_not_index = vec(readdlm("julia_merger_result/subsidy_used_or_not_index_shared.txt",',',Float64))
res_theta_solved = copy(res_theta)
for i = 1:size(res_theta)[1], j = 1:size(res_theta)[2]
	if res_theta[i,j] == 0.0
		res_theta_solved[i,j] = copy(NaN)
	end
end

res_correct_ineq./num_all_ineq
Plots.histogram(res_theta_solved[subsidy_used_or_not_index.<1,1],
                bins = 60,
                ylabel = "count",
                xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated β",
				title = "subsidy (shared)",
				xlim = [-20,20],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[subsidy_used_or_not_index.>0,1],
                bins = 60,
                label = "estimated β (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.β_0], label="true β", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_beta_shared")

Plots.histogram(res_theta_solved[:,2],
                bins = 60,
                ylabel = "count",
				xlabel = "δ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
				label = "estimated δ",
				title = "subsidy (shared)",
				xlim = [-20,20],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[subsidy_used_or_not_index.>0,2],
                bins = 60,
                label = "estimated δ (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.δ_0], label="true δ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_delta_shared")

Plots.histogram(res_theta_solved[:,3],
                bins = 60,
                ylabel = "count",
				xlabel = "γ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated γ",
				title = "subsidy (shared)",
				xlim = [-20,20],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[subsidy_used_or_not_index.>0,3],
                bins = 60,
                label = "estimated γ (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_gamma_shared")

median_bias_beta_with_subsidy_shared = nanmedian(res_theta_solved[subsidy_used_or_not_index.<1,1].-m.β_0)
median_bias_beta_without_subsidy_shared = nanmedian(res_theta_solved[subsidy_used_or_not_index.>0,1].-m.β_0)
median_bias_gamma_with_subsidy_shared = nanmedian(res_theta_solved[subsidy_used_or_not_index.<1,3].-m.γ_0)
median_bias_gamma_without_subsidy_shared = nanmedian(res_theta_solved[subsidy_used_or_not_index.>0,3].-m.γ_0)
RMSE_beta_with_subsidy_shared = RMSE(res_theta_solved[subsidy_used_or_not_index.<1,1],m.β_0)
RMSE_beta_without_subsidy_shared = RMSE(res_theta_solved[subsidy_used_or_not_index.>0,1],m.β_0)
RMSE_gamma_with_subsidy_shared = RMSE(res_theta_solved[subsidy_used_or_not_index.<1,3],m.γ_0)
RMSE_gamma_without_subsidy_shared = RMSE(res_theta_solved[subsidy_used_or_not_index.>0,3],m.γ_0)

#----------------------------------------------#
# create monte carlo table
#----------------------------------------------#
LaTeXTabulars.latex_tabular("julia_merger_table/monte_carlo_results_bias_RMSE.tex",
			  Tabular("@{\\extracolsep{5pt}}ccccc"),
			  [Rule(:top),
			   ["Num of firms", "subsidy type", "parameter",
			   "median error", "RMSE"],
			   Rule(:mid),
			   ["$(m.N)", "to buyer", L"\beta",
			   round(median_bias_beta_with_subsidy_to_buyer,digits=2),
			   round(RMSE_beta_with_subsidy_to_buyer,digits=2)],
			   ["$(m.N)", "to buyer", L"\gamma",
			   round(median_bias_gamma_with_subsidy_to_buyer,digits=2),
			   round(RMSE_gamma_with_subsidy_to_buyer,digits=2)],
			   ["$(m.N)", "shared", L"\beta",
 			   round(median_bias_beta_with_subsidy_shared,digits=2),
 			   round(RMSE_beta_with_subsidy_shared,digits=2)],
 			   ["$(m.N)", "shared", L"\gamma",
 			   round(median_bias_gamma_with_subsidy_shared,digits=2),
 			   round(RMSE_gamma_with_subsidy_shared,digits=2)],
			   Rule(:bottom)]
			   )
#-------------------------------------#
# identification power with no merger
#-------------------------------------#

temp_threshold_tonnage = 1
temp_subsidy_amount = 1
temp_subsidy_type = "shared"
@time m = carrier_mod(randomseed = 3,#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 15.0,# coefficient of additional merger cost
					 ton_dim= 2)
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = temp_subsidy_type)
@show obsd.tarid
temp_simulate_num = 1000
#----------------------------------------------#
# plot subsidy (to buyer)
#----------------------------------------------#
you_want_run = "run"
if you_want_run == "run"
    @time res_theta, ~, ~, ~, ~, ~, ~, all_unmatched_or_not_index = maxscore_mc_temp(8,
	                     num_its = temp_simulate_num,
						 threshold_tonnage = 1,
						 subsidy_amount = 1,
						 subsidy_type = temp_subsidy_type,
						 β_0_temp = 0.0,
						 δ_0_temp = 1.0,
						 γ_0_temp = 15.0,
						 search_range_temp = (-10.0,40.0))
	open("julia_merger_result/res_theta_no_merger_$temp_subsidy_type.txt", "w") do io
		DelimitedFiles.writedlm(io, res_theta,",")
	end
	open("julia_merger_result/all_unmatched_or_not_index_no_merger_$temp_subsidy_type.txt", "w") do io
		DelimitedFiles.writedlm(io, all_unmatched_or_not_index,",")
	end
end


res_theta = readdlm("julia_merger_result/res_theta_no_merger_$temp_subsidy_type.txt",',',Float64)
all_unmatched_or_not_index = vec(readdlm("julia_merger_result/all_unmatched_or_not_index_no_merger_$temp_subsidy_type.txt",',',Float64))
res_theta_solved = copy(res_theta)
for i = 1:size(res_theta)[1], j = 1:size(res_theta)[2]
	if res_theta[i,j] == 0.0
		res_theta_solved[i,j] = copy(NaN)
	end
end

#res_correct_ineq./num_all_ineq
Plots.histogram(res_theta_solved[all_unmatched_or_not_index.<1,1],
                bins = 60,
                ylabel = "count",
                xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated β",
				title = "subsidy ($temp_subsidy_type)",
				xlim = [-10,40],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_or_not_index.>0,1],
                bins = 60,
                label = "estimated β (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.β_0], label="true β", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_beta_no_merger_$temp_subsidy_type")

Plots.histogram(res_theta_solved[all_unmatched_or_not_index.<1,3],
                bins = 60,
                ylabel = "count",
				xlabel = "γ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated γ",
				title = "subsidy ($temp_subsidy_type)",
				xlim = [-10,40],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_or_not_index.>0,3],
                bins = 60,
                label = "estimated γ (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_gamma_no_merger_$temp_subsidy_type")

#-----------------------#
# plot to_buyer
#-----------------------#
temp_subsidy_type = "to_buyer"

you_want_run = "run"
if you_want_run == "run"
    @time res_theta, ~, ~, ~, ~, ~, ~, all_unmatched_or_not_index = maxscore_mc_temp(8,
	                     num_its = temp_simulate_num,
						 threshold_tonnage = 1,
						 subsidy_amount = 1,
						 subsidy_type = temp_subsidy_type,
						 β_0_temp = 0.0,
						 δ_0_temp = 1.0,
						 γ_0_temp = 15.0,
						 search_range_temp = (-10.0,40.0))
	open("julia_merger_result/res_theta_no_merger_$temp_subsidy_type.txt", "w") do io
		DelimitedFiles.writedlm(io, res_theta,",")
	end
	open("julia_merger_result/all_unmatched_or_not_index_no_merger_$temp_subsidy_type.txt", "w") do io
		DelimitedFiles.writedlm(io, all_unmatched_or_not_index,",")
	end
end


res_theta = readdlm("julia_merger_result/res_theta_no_merger_$temp_subsidy_type.txt",',',Float64)
all_unmatched_or_not_index = vec(readdlm("julia_merger_result/all_unmatched_or_not_index_no_merger_$temp_subsidy_type.txt",',',Float64))
res_theta_solved = copy(res_theta)
for i = 1:size(res_theta)[1], j = 1:size(res_theta)[2]
	if res_theta[i,j] == 0.0
		res_theta_solved[i,j] = copy(NaN)
	end
end

#res_correct_ineq./num_all_ineq
Plots.histogram(res_theta_solved[all_unmatched_or_not_index.<1,1],
                bins = 60,
                ylabel = "count",
                xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated β",
				title = "subsidy ($temp_subsidy_type)",
				xlim = [-10,40],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_or_not_index.>0,1],
                bins = 60,
                label = "estimated β (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.β_0], label="true β", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_beta_no_merger_$temp_subsidy_type")

Plots.histogram(res_theta_solved[all_unmatched_or_not_index.<1,3],
                bins = 60,
                ylabel = "count",
				xlabel = "γ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated γ",
				title = "subsidy ($temp_subsidy_type)",
				xlim = [-10,40],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_or_not_index.>0,3],
                bins = 60,
                label = "estimated γ (subsidy not used)",
				alpha = 0.3)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_gamma_no_merger_$temp_subsidy_type")

#--------------------------#
# plot score function (no-merger market)
#--------------------------#
@time m = carrier_mod(randomseed = 4,#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 15.0,# coefficient of additional merger cost
					 ton_dim= 2)
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = temp_subsidy_type)
m.Bundle[obsd.tarid]
obsd.total_size_target + obsd.total_size_buyer
total_unmatched = sum(obsd.tarid .== 1)
total_num_group = length(unique(obsd.tarid[obsd.tarid .!= 1]))

domain1 = [-20:2.0:30;]
domain2 = [-20:2.0:30;]
res_contor = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor[i,j], ~ = score_b(m, [domain1[i], m.δ_0, domain2[j]],
	                             obsd, m.N, temp_threshold_tonnage,
								 temp_subsidy_amount,
								 subsidy_type = temp_subsidy_type)
end
Plots.contour(domain1, domain2, res_contor', fill = true,
              ylabel = "γ (merger cost)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)], $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
			  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
Plots.hline!([m.γ_0], label="true γ", linestyle = :dash)
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_gamma_no_merger_$(temp_subsidy_type)_subsidy")


#--------------------------#
# duopoly and oligopoly market
#--------------------------#

N_list = [2:1:7;]
seed_list = [4,4,4,5,4,6] # to show illustrative bounds
@time for ii = 1:length(N_list)
	@show N_iter = N_list[ii]
	@time m = carrier_mod(randomseed = seed_list[ii],#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = N_iter,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)
	utility, obsd = gen_data(m,
	                threshold_tonnage =　temp_threshold_tonnage,
	                subsidy_amount = temp_subsidy_amount,
					subsidy_type = temp_subsidy_type)
	@show m.Bundle[obsd.tarid]
	obsd.total_size_target + obsd.total_size_buyer
	total_unmatched = sum(obsd.tarid .== 1)
	total_num_group = length(unique(obsd.tarid[obsd.tarid .!= 1]))
	domain1 = [-20:2.0:30;]
	domain2 = [-20:2.0:30;]
	res_contor = zeros(length(domain1),length(domain2))
	for i in 1:length(domain1),j in 1:length(domain2)
		res_contor[i,j], ~ = score_b(m, [domain1[i], m.δ_0, domain2[j]],
		                             obsd, m.N, temp_threshold_tonnage,
									 temp_subsidy_amount,
									 subsidy_type = temp_subsidy_type)
	end
	Plots.contour(domain1, domain2, res_contor', fill = true,
	              ylabel = "γ (merger cost)",
				  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)], $(total_unmatched)/$(m.N) unmatched, $(total_num_group) group(s)",
				  title = "Objective value: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
	Plots.hline!([m.γ_0], label="true γ", linestyle = :dash)
	Plots.vline!([m.β_0], label="true β", linestyle = :dash)
	N_iter = Int(N_iter)
	savefig("julia_merger_figure/contour_beta_gamma_N_$(N_iter)_$(temp_subsidy_type)_subsidy")
end
