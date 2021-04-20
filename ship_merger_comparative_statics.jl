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
using JLD # restore for multidimensional array
include(joinpath(dirname(@__FILE__),"functions.jl"))
using .functions
temp_randomseed = 1
#---------------------#
# comparative statics
#---------------------#
@time m = carrier_mod(randomseed = 4,#1
					 ton_dist = Distributions.LogNormal(2.0,1.0),
					 N = 8,
					 β_0 = 0.0,# coefficient of covariates
					 δ_0 = 1.0,# coefficient of subsidy
					 γ_0 = 0.0,# coefficient of additional merger cost
					 ton_dim= 2)
#-----------------------------------#
# (1) calculate the number of firms
#-----------------------------------#
temp_threshold_tonnage = 1
temp_subsidy_amount = 1
temp_subsidy_type = "to_buyer"
function gen_plot_number_of_firms(m,
	                             temp_threshold_tonnage,
								 temp_subsidy_amount,
								 temp_subsidy_type;
								 seed_max = 100)
	utility, obsd = gen_data(m,
					threshold_tonnage = temp_threshold_tonnage,
					subsidy_amount = temp_subsidy_amount,
					subsidy_type = temp_subsidy_type)
	domain1 = [0:1.0:10;]
	domain2 = [0:1.0:5;]
	num_of_group = ones(length(domain1), length(domain2), seed_max).*100
	num_of_unmatched = ones(length(domain1), length(domain2), seed_max).*100
	num_of_group_and_unmatched = ones(length(domain1), length(domain2), seed_max).*100
	num_of_target_id = ones(length(domain1), length(domain2), seed_max, m.N).*100
	for i in 1:length(domain1), j in 1:length(domain2), k in 1:seed_max
		iter = domain1[i]
		iter2 = domain2[j]
		@time m = carrier_mod(randomseed = k,#1
		                     ton_dist = Distributions.LogNormal(2.0,1.0),
							 N = 8,
							 β_0 = 0.0,# coefficient of covariates
							 δ_0 = iter2,# coefficient of subsidy
							 γ_0 = iter,# coefficient of additional merger cost
							 ton_dim= 2)
		utility, obsd = gen_data(m,
		                threshold_tonnage = temp_threshold_tonnage,
		                subsidy_amount = temp_subsidy_amount,
						subsidy_type = temp_subsidy_type)
		if length(obsd.tarid) == m.N
			# integer matching
			global buyer_index = zeros(0)
			for ii = 1:m.N
				if IS_Sell(m, ii, obsd.tarid[ii]) != 1 && obsd.tarid[ii] !== 1
					global buyer_index = vcat(buyer_index, ii)
				end
			end
			num_of_target_id[i,j,k,:] = obsd.tarid
			num_of_group[i,j,k] = length(buyer_index)
			num_of_unmatched[i,j,k] = length(obsd.tarid[obsd.tarid .== 1])
			num_of_group_and_unmatched[i,j,k] = num_of_group[i,j,k] + num_of_unmatched[i,j,k]
		end
	end
	mean_num_of_group = zeros(length(domain1),length(domain2))
	mean_num_of_unmatched = zeros(length(domain1),length(domain2))
	mean_num_of_group_and_unmatched = zeros(length(domain1),length(domain2))
	for i in 1:length(domain1), j in 1:length(domain2)
		mean_num_of_group[i,j] = mean(num_of_group[i,j,:][num_of_group[i,j,:].!=100])
		mean_num_of_unmatched[i,j] = mean(num_of_unmatched[i,j,:][num_of_unmatched[i,j,:].!=100])
		mean_num_of_group_and_unmatched[i,j] = mean_num_of_group[i,j] + mean_num_of_unmatched[i,j]
	end

	Plots.plot(domain1, mean_num_of_group[:,1],
	           label="Num of groups [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[1]), γ]",
	           xlabel="γ",
			   ylabel="number of groups",
			   legend=:topright,
			   linestyle = :dash,
			   alpha = 0.7,
			   markershape = :auto,
	           title = "Num of groups: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")
	for ii = 2:length(mean_num_of_group[1,:])
		Plots.plot!(domain1, mean_num_of_group[:,ii],
		           linestyle = :dash,
		           alpha = 0.7,
		           markershape = :auto,
		           label="Num of groups [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[ii]), γ]")
	end
	Plots.plot!()
	savefig("julia_merger_figure/plot_num_of_groups_$(temp_subsidy_type)_subsidy")

	Plots.plot(domain1, mean_num_of_unmatched[:,1],
	           label="Num of unmatched [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[1]), γ]",
	           xlabel="γ",
			   ylabel="number of unmatched",
			   legend=:bottomright,
			   linestyle = :dash,
			   alpha = 0.7,
			   ylim = [0,m.N],
			   markershape = :auto,
	           title = "Num of unmatched: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")

	for ii = 2:length(mean_num_of_unmatched[1,:])
		Plots.plot!(domain1, mean_num_of_unmatched[:,ii],
		           linestyle = :dash,
				   alpha = 0.7,
				   markershape = :auto,
		           label="Num of unmatched [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[ii]), γ]")
	end
	Plots.hline!([m.N], label="",
	             linestyle = :dash, color = :black)
	Plots.plot!()
	savefig("julia_merger_figure/plot_num_of_unmatched_$(temp_subsidy_type)_subsidy")

	Plots.plot(domain1, mean_num_of_group_and_unmatched[:,1],
	           label="Num of post-merger firms [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[1]), γ]",
	           xlabel="γ",
			   ylabel="number of post-merger firms",
			   legend=:bottomright,
			   linestyle = :dash,
			   alpha = 0.7,
			   ylim = [0,m.N],
			   markershape = :auto,
	           title = "Num of post-merger firms: (threshold, amount, type)=($(temp_threshold_tonnage),$(temp_subsidy_amount),$(temp_subsidy_type))")

	for ii = 2:length(mean_num_of_group_and_unmatched[1,:])
		Plots.plot!(domain1, mean_num_of_group_and_unmatched[:,ii],
		           linestyle = :dash,
				   alpha = 0.7,
				   markershape = :auto,
		           label="Num of post-merger firms [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[ii]), γ]")
	end
	Plots.hline!([m.N], label="",
	             linestyle = :dash, color = :black)
	Plots.plot!()
	savefig("julia_merger_figure/plot_num_of_post_merger_firms_$(temp_subsidy_type)_subsidy")
	return num_of_target_id
end
temp_subsidy_type = "to_buyer"
@time num_of_target_id = gen_plot_number_of_firms(m,
                         temp_threshold_tonnage,
						 temp_subsidy_amount,
						 temp_subsidy_type;
						 seed_max = 50)
JLD.save("julia_merger_result/num_of_target_id_$(temp_subsidy_type)_subsidy.jld",
         "data", num_of_target_id)
num_of_target_id = JLD.load("julia_merger_result/num_of_target_id_$(temp_subsidy_type)_subsidy.jld")["data"]

#-------------------------#
# shared subsidy case
#-------------------------#
temp_threshold_tonnage = 1
temp_subsidy_amount = 1
temp_subsidy_type = "shared"
@time num_of_target_id = gen_plot_number_of_firms(m,
                         temp_threshold_tonnage,
						 temp_subsidy_amount,
						 temp_subsidy_type;
						 seed_max = 50)
JLD.save("julia_merger_result/num_of_target_id_$(temp_subsidy_type)_subsidy.jld",
         "data", num_of_target_id)
num_of_target_id = JLD.load("julia_merger_result/num_of_target_id_$(temp_subsidy_type)_subsidy.jld")["data"]
#----------------------------------------#
# (2) calculate the merger configuration
#----------------------------------------#


domain1 = vcat([0:1.0:10;],20)
model_list = zeros(length(domain1),length(domain2), m.N)
k = 8
for i in 1:length(domain1), j in 1:length(domain2)
	iter1 = domain1[i]
	iter2 = domain2[j]
	@time m = carrier_mod(randomseed = k,#1
						 ton_dist = Distributions.LogNormal(2.0,1.0),
						 N = 8,
						 β_0 = 0.0,# coefficient of covariates
						 δ_0 = iter2,# coefficient of subsidy
						 γ_0 = iter1,# coefficient of additional merger cost
						 ton_dim= 2)
	utility, obsd = gen_data(m,
					threshold_tonnage = temp_threshold_tonnage,
					subsidy_amount = temp_subsidy_amount,
					subsidy_type = temp_subsidy_type)
	if length(obsd.tarid) == m.N
		# integer equilibrium
		global total_tonnage_market = sum(obsd.total_size_buyer)
		buyer_index = zeros(0)
		unmatched_index = zeros(0)
		for ii = 1:m.N
			if IS_Sell(m, ii, obsd.tarid[ii]) != 1 && obsd.tarid[ii] !== 1
				buyer_index = vcat(buyer_index, ii)
			end
			if obsd.tarid[ii] == 1
				unmatched_index = vcat(unmatched_index, ii)
			end
		end
		buyer_index = Int.(buyer_index)
		unmatched_index = Int.(unmatched_index)
		seller_length = m.N -length(buyer_index) - length(unmatched_index)
		model_list[i,j,:] = vcat(obsd.total_size_target[buyer_index] + obsd.total_size_buyer[buyer_index],
								 obsd.total_size_buyer[unmatched_index],
								 zeros(seller_length))
	else
		model_list[i,j,:] = vcat(zeros(m.N))
	end
end
obsd.total_size_buyer
using StatsPlots
model_list_ordered = zeros(length(domain1),length(domain2), m.N)
for i in 1:length(domain1), j in 1:length(domain2)
	model_list_ordered[i,j,:] = sort(model_list[i,j,:],rev=true)
end
for jj = 1:length(domain2)
	StatsPlots.groupedbar(
	    sortslices(model_list_ordered[:,jj,:],dims=2,rev=true),# sort by group-size
        bar_position = :stack,
        bar_width=0.7,
		title="Merger configuration: [β₀,δ₀,γ₀] = [$(m.β_0), $(domain2[jj]), γ],$(temp_subsidy_type) subsidy",
		xlabel="γ",
		ylabel="Total firm size",
		legend = false,
        xticks=(1:length(domain1), domain1)
		)
	Plots.hline!([total_tonnage_market], label="",
	             linestyle = :dash, color = :black)
	temp = Int(domain2[jj])
	savefig("julia_merger_figure/stacked_plot_merger_configuration_delta_$(temp)_$(temp_subsidy_type)_subsidy")
end







# end
