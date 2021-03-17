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
#-----------------------#
# setup
#-----------------------#
@time m = carrier_mod(randomseed = 8,
					 N = 8,
					 β_0 = 3.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 10.0,# coefficient of additional merger cost
					 ton_dim= 2)
#utility, obsd = gen_data(m, threshold_tonnage = 0, subsidy_amount = 0)
#utility, obsd = gen_data(m, threshold_tonnage =　70, subsidy_amount = 0)
#utility, obsd = gen_data(m, threshold_tonnage =　70, subsidy_amount = 50)
temp_threshold_tonnage = 70
temp_subsidy_amount = 1
temp_subsidy_type = "shared"
#utility, obsd = gen_data(m, threshold_tonnage =　temp_threshold_tonnage, subsidy_amount = 0, subsidy_type = temp_subsidy_type)
utility, obsd = gen_data(m, threshold_tonnage =　temp_threshold_tonnage, subsidy_amount = temp_subsidy_amount, subsidy_type = temp_subsidy_type)
~,~,temp = score_b(m, [m.β_0, m.δ_0, m.γ_0], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = temp_subsidy_type)
#~,~,temp = score_b(m, [m.β_0, m.δ_0, 50], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = temp_subsidy_type)

domain = [-10:1.0:50;]
res_β_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_β_0[i], ~ = score_b(m, [iter, m.δ_0, m.γ_0], obsd,
	                        m.N, temp_threshold_tonnage,
							temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_β_0, label="score", xlabel="beta",title = "Objective value: [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]")
Plots.vline!([m.β_0], label="true", linestyle = :dot)

domain = [-10:1.0:10;]
res_δ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_δ_0[i], ~ = score_b(m, [m.β_0, iter, m.γ_0], obsd,
	                        m.N, temp_threshold_tonnage,
							temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_δ_0, label="score", xlabel="delta(coefficient of subsidy)",title = "Objective value: [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]")
Plots.vline!([m.δ_0], label="true", linestyle = :dot)

domain = [-10:1.0:100;]
res_γ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_γ_0[i], ~ = score_b(m, [m.β_0, m.δ_0, iter], obsd,
	                        m.N, temp_threshold_tonnage,
							temp_subsidy_amount,
							subsidy_type = temp_subsidy_type)
end
Plots.plot(domain, res_γ_0, label="score", xlabel="gamma(coefficient of additional merger cost)",title = "Objective value: [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]")
Plots.vline!([m.γ_0], label="true", linestyle = :dot)

#--------------------------#
# Identification plots
#--------------------------#
@time m = carrier_mod(randomseed = 3,#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 4.0,# coefficient of covariates
					 δ_0 = 2.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)
@show m.ton
temp_threshold_tonnage = 100
temp_subsidy_amount = 10
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = "shared")
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
domain = [-10:1.0:50;]
res_β_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_β_0[i], ~ = score_b(m, [iter, m.δ_0, m.γ_0], obsd,
	                m.N, temp_threshold_tonnage,
					temp_subsidy_amount,
					subsidy_type = "to_buyer")
end
Plots.plot(domain, res_β_0, label="score",
           xlabel="β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
           title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.vline!([m.β_0], label="true", linestyle = :dot)
savefig("julia_merger_figure/plot_beta")

domain = [-10:1.0:50;]
res_δ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_δ_0[i], ~ = score_b(m, [m.β_0, iter, m.γ_0], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.plot(domain, res_δ_0, label="score",
           xlabel="δ (coefficient of subsidy) where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
           title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount)) ")
Plots.vline!([m.δ_0], label="true", linestyle = :dot)
savefig("julia_merger_figure/plot_delta")

domain = [-10:1.0:50;]
res_γ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_γ_0[i], ~ = score_b(m, [m.β_0, m.δ_0, iter], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.plot(domain, res_γ_0, label="score",
           xlabel = "γ (coefficient of additional merger cost) where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
           title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount)) ")
Plots.vline!([m.γ_0], label="true", linestyle = :dot)
savefig("julia_merger_figure/plot_gamma")

domain1 = [-10:1.0:20;]
domain2 = [-10:1.0:20;]
res_contor = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor[i,j], ~ = score_b(m, [domain1[i], m.δ_0, domain2[j]], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.contour(domain1, domain2, res_contor, fill = true,
              ylabel = "γ (coefficient of additional merger cost)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
			  title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.hline!([m.γ_0], label="true γ", linestyle = :dash)
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_gamma")

domain1 = [-10:1.0:20;]
domain2 = [-10:1.0:20;]
res_contor2 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor2[i,j], ~ = score_b(m, [m.β_0, domain1[i], domain2[j]], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.contour(domain1, domain2, res_contor2, fill = true,
              ylabel = "γ (coefficient of additional merger cost)",
			  xlabel = "δ (coefficient of subsidy) where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
			  title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash)
savefig("julia_merger_figure/contour_gamma_delta")

domain1 = [-10:1.0:20;]
domain2 = [-10:1.0:20;]
res_contor3 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor3[i,j], ~ = score_b(m, [domain1[i], domain2[j], m.γ_0], obsd, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.contour(domain1, domain2, res_contor3, fill = true,
              ylabel = "δ (coefficient of subsidy)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
			  title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_delta")

#--------------------------#
# (2) s_to_buyer
#--------------------------#

domain = [-10:1.0:50;]
res_β_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_β_0[i], ~ = score_b(m, [iter, m.δ_0, m.γ_0], obsd_to_buyer,
	                m.N, temp_threshold_tonnage,
					temp_subsidy_amount,
					subsidy_type = "to_buyer")
end
Plots.plot(domain, res_β_0, label="score",
           xlabel="β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
           title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.vline!([m.β_0], label="true", linestyle = :dot)
savefig("julia_merger_figure/plot_beta_to_buyer")

domain = [-10:1.0:50;]
res_δ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_δ_0[i], ~ = score_b(m, [m.β_0, iter, m.γ_0], obsd_to_buyer,
	 m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.plot(domain, res_δ_0, label="score",
           xlabel="δ (coefficient of subsidy) where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
           title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount)) ")
Plots.vline!([m.δ_0], label="true", linestyle = :dot)
savefig("julia_merger_figure/plot_delta_to_buyer")

domain = [-10:1.0:50;]
res_γ_0 = zeros(length(domain))
for i in 1:length(domain)
	iter = domain[i]
	res_γ_0[i], ~ = score_b(m, [m.β_0, m.δ_0, iter], obsd_to_buyer,
	 m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.plot(domain, res_γ_0, label="score",
           xlabel = "γ (coefficient of additional merger cost) where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
           title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount)) ")
Plots.vline!([m.γ_0], label="true", linestyle = :dot)
savefig("julia_merger_figure/plot_gamma_to_buyer")

domain1 = [-10:1.0:20;]
domain2 = [-10:1.0:20;]
res_contor = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor[i,j], ~ = score_b(m, [domain1[i], m.δ_0, domain2[j]],
	                             obsd_to_buyer, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.contour(domain1, domain2, res_contor, fill = true,
              ylabel = "γ (coefficient of additional merger cost)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
			  title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.hline!([m.γ_0], label="true γ", linestyle = :dash)
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_gamma_to_buyer")

domain1 = [-10:1.0:20;]
domain2 = [-10:1.0:20;]
res_contor2 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor2[i,j], ~ = score_b(m, [m.β_0, domain1[i], domain2[j]],
	                              obsd_to_buyer, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.contour(domain1, domain2, res_contor2, fill = true,
              ylabel = "γ (coefficient of additional merger cost)",
			  xlabel = "δ (coefficient of subsidy) where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
			  title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash)
savefig("julia_merger_figure/contour_gamma_delta_to_buyer")

domain1 = [-10:1.0:20;]
domain2 = [-10:1.0:20;]
res_contor3 = zeros(length(domain1),length(domain2))
for i in 1:length(domain1),j in 1:length(domain2)
	res_contor3[i,j], ~ = score_b(m, [domain1[i], domain2[j], m.γ_0],
	                              obsd_to_buyer, m.N, temp_threshold_tonnage, temp_subsidy_amount, subsidy_type = "to_buyer")
end
Plots.contour(domain1, domain2, res_contor3, fill = true,
              ylabel = "δ (coefficient of subsidy)",
			  xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
			  title = "Objective value: (threshold, amount)=($(temp_threshold_tonnage),$(temp_subsidy_amount))")
Plots.vline!([m.β_0], label="true β", linestyle = :dash)
Plots.hline!([m.δ_0], label="true δ", linestyle = :dash)
savefig("julia_merger_figure/contour_beta_delta_to_buyer")




#--------------------------#
# (3)
#--------------------------#
function maxscore_mc_temp(num_agents::Int64;
                     num_its::Int64 = 10,
                     threshold_tonnage = 100,
					 subsidy_amount = 10,
					 subsidy_type = "to_buyer")
    start = Dates.unix2datetime(time())
	param_dim = 3
    myests = zeros(num_its, param_dim)
	mynum_correct_ineq = zeros(num_its, 1)
	truenum_correct_ineq = zeros(num_its, 1)
	num_all_ineq = zeros(num_its)
	solved_case_index = zeros(Int64, 0)
	all_unmatched_index = zeros(num_its)
    for i = 1:num_its
	   @show i
       println("Create obsdat for iteration $i \n" )
	   @time m = carrier_mod(randomseed = i,
							N = num_agents,
							β_0 = 4.0,# coefficient of covariates
	  					    δ_0 = 2.0,# coefficient of subsidy
	  					    γ_0 = 4.0,# coefficient of additional merger cost
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
		   if sum(obsd.tarid .!=1)==0
			   println("simulated data has only unmatched firms")
			   all_unmatched_index[i] = 1
		   end
		   solved_case_index = vcat(solved_case_index, i)
		   function score_bthis(beta::Vector{Float64}, obsdat::DataFrame)
			   score, ~, ~ = score_b(m, beta, obsdat, num_agents, threshold_tonnage, subsidy_amount, subsidy_type = subsidy_type)
			   res = -1.0*score + 100000.0 # need to be Float64 for bboptimize
			   println("score:$(res) \n" )
			   return res
		   end
		   m_res = BlackBoxOptim.bboptimize(beta -> score_bthis(beta, obsd);
		                                    SearchRange = (-10.0, 10.0),
		                                    NumDimensions = param_dim, Method = :de_rand_1_bin,  MaxSteps=50)
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
   all_unmatched_index
   #@show myests
   return myests, res_abs_mean_err, res_sqrt, mynum_correct_ineq, truenum_correct_ineq, num_all_ineq, all_unmatched_index
end


@time m = carrier_mod(randomseed = 1,
					 N = 8,
					 β_0 = 4.0,# coefficient of covariates
					 δ_0 = 2.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)

temp_simulate_num = 1000
#----------------------------------------------#
# plot subsidy (to buyer)
#----------------------------------------------#
you_want_run = "not run"
if you_want_run == "run"
    @time res_theta, res_abs_mean_err, res_sqrt, res_correct_ineq, res_truenum_correct_ineq, num_all_ineq, all_unmatched_index = maxscore_mc_temp(8, num_its = temp_simulate_num,
                                                                                                                      threshold_tonnage = 100,
																													  subsidy_amount = 10,
																													  subsidy_type = "to_buyer")
	open("julia_merger_result/res_theta_to_buyer.txt", "w") do io
		DelimitedFiles.writedlm(io, res_theta,",")
	end
	open("julia_merger_result/all_unmatched_index_to_buyer.txt", "w") do io
		DelimitedFiles.writedlm(io, all_unmatched_index,",")
	end
end
res_theta = readdlm("julia_merger_result/res_theta_to_buyer.txt",',',Float64)
all_unmatched_index = vec(readdlm("julia_merger_result/all_unmatched_index_to_buyer.txt",',',Float64))
#res_theta_solved = res_theta[res_theta.!=0]
res_theta_solved = copy(res_theta)
for i = 1:size(res_theta)[1], j = 1:size(res_theta)[2]
	if res_theta[i,j] == 0.0
		res_theta_solved[i,j] = copy(NaN)
	end
end
Plots.histogram(res_theta_solved[all_unmatched_index.<1,1], bins=40,
                ylabel = "count",
                xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated β",
				title = "subsidy (to buyer)",
				xlim =[-10,10],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_index.>0,1], bins=40,
                label = "estimated β (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.β_0], label="true β", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_beta_to_buyer")

Plots.histogram(res_theta_solved[all_unmatched_index.<1,2], bins=40,
                ylabel = "count",
				xlabel = "δ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
				label = "estimated δ",
				title = "subsidy (to buyer)",
				xlim =[-10,10],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_index.>0,2], bins=40,
                label = "estimated δ (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.δ_0], label="true δ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_delta_to_buyer")

Plots.histogram(res_theta_solved[all_unmatched_index.<1,3], bins=40,
                ylabel = "count",
				xlabel = "γ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated γ",
				title = "subsidy (to buyer)",
				xlim =[-10,10],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_index.>0,3], bins=40,
                label = "estimated γ (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_gamma_to_buyer")

#----------------------------------------------#
# plot subsidy (shared)
#----------------------------------------------#
you_want_run = "not run"
if you_want_run == "run"
	@time res_theta, res_abs_mean_err, res_sqrt, res_correct_ineq, res_truenum_correct_ineq, num_all_ineq, all_unmatched_index = maxscore_mc_temp(8, num_its = temp_simulate_num,
	                                                                                                                      threshold_tonnage = 100,
																														  subsidy_amount = 10,
																														  subsidy_type = "shared")
	open("julia_merger_result/res_theta_shared.txt", "w") do io
		DelimitedFiles.writedlm(io, res_theta,",")
	end
	open("julia_merger_result/all_unmatched_index_shared.txt", "w") do io
		DelimitedFiles.writedlm(io, all_unmatched_index,",")
	end
end

res_theta = readdlm("julia_merger_result/res_theta_shared.txt",',',Float64)
all_unmatched_index = vec(readdlm("julia_merger_result/all_unmatched_index_shared.txt",',',Float64))
res_theta_solved = copy(res_theta)
for i = 1:size(res_theta)[1], j = 1:size(res_theta)[2]
	if res_theta[i,j] == 0.0
		res_theta_solved[i,j] = copy(NaN)
	end
end

res_correct_ineq./num_all_ineq
Plots.histogram(res_theta_solved[all_unmatched_index.<1,1], bins=40,
                ylabel = "count",
                xlabel = "β where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated β",
				title = "subsidy (shared)",
				xlim = [-10,10],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_index.>0,1], bins=40,
                label = "estimated β (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.β_0], label="true β", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_beta_shared")

Plots.histogram(res_theta_solved[:,2], bins=40,
                ylabel = "count",
				xlabel = "δ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
				label = "estimated δ",
				title = "subsidy (shared)",
				xlim =[-10,10],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_index.>0,2], bins=40,
                label = "estimated δ (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.δ_0], label="true δ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_delta_shared")

Plots.histogram(res_theta_solved[:,3], bins=40,
                ylabel = "count",
				xlabel = "γ where [β₀,δ₀,γ₀] = [$(m.β_0), $(m.δ_0), $(m.γ_0)]",
                label = "estimated γ",
				title = "subsidy (shared)",
				xlim =[-10,10],
				alpha = 0.3)
Plots.histogram!(res_theta_solved[all_unmatched_index.>0,3], bins=40,
                label = "estimated γ (all unmatched)",
				alpha = 0.3)
Plots.vline!([m.γ_0], label="true γ", linestyle = :dash,
             color = "black")
savefig("julia_merger_figure/histogram_gamma_shared")


#------------------------#
# comparative statistics
#------------------------#
@time m = carrier_mod(randomseed = 3,#1
                     ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 4.0,# coefficient of covariates
					 δ_0 = 2.0,# coefficient of subsidy
					 γ_0 = 4.0,# coefficient of additional merger cost
					 ton_dim= 2)
temp_threshold_tonnage = 100
temp_subsidy_amount = 10
utility, obsd = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = "shared")
utility_to_buyer, obsd_to_buyer = gen_data(m,
                threshold_tonnage =　temp_threshold_tonnage,
                subsidy_amount = temp_subsidy_amount,
				subsidy_type = "to_buyer")
@show obsd.tarid
@show obsd_to_buyer.tarid
m.Bundle[obsd_to_buyer.tarid]
