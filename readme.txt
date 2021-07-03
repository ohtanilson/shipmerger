# Estimating Endogenous Coalitional Mergers: Assortativeness of Size and Techonogical Specialization.
# written by Suguru Otani, 2021.06.30

PATH must be on shipmerger/.
The following folders contain all R implementations and raw data.
- shipmerger/main: This contains execution files.
- shipmerger/input: This contains a raw dataset.
- shipmerger/output: This contains cleaned dataset for Julia.
- shipmerger/figuretable: This contains figures and tables shown in my paper.

Since computing all processes of Julia on a single laptop may take a few weeks,
I locate the following folders restoring current Julia result files on ../shipmerger/.
- shipmerger/julia_merger_result : This contains all intermediate files.
- shipmerger/julia_merger_table: This contains all tables converting intermediate files into latex style
- shipmerger/julia_merger_figure: This contains all figures.

The following files restore (and update current files if you run a code) the results on these folders. 
You can replicate all results by following processes sequentially.

#---------------------------------------------------#
# Data Cleaning and Reduced Form Analysis (Rproject)
#---------------------------------------------------#
- shipmerger/main/01_make_figuretable.R
  This cleans raw data file (210627firmdata_processed.csv) in shipmerger/input.
  Then, this produces cleaned data file (data_for_maximum_rank_estimation.csv) for Julia in shipmerger/output.
  This also generates related Figures and Tables in my paper on shipmerger/figuretable.
  It takes just a minute.
- shipmerger/main/02_make_sankey_diagram
  This generates sankey diagram in Figure 9. 
  Input data is manually assigned based on results on ship_merger_counterfactual.jl.


#-----------------------------------------#
# Matching Maximum Rank Estimation (Julia)
#-----------------------------------------#
- shipmerger/ship_merger_score_estimation.jl
  This computes 8 specifications with bootstrap 95% CI for
    (1). multivariate matching maximum rank estimation in Table 10
         (22 hours for 200 point-estimates, 21 hours for CI bootstrapping 200)
    (2). matching maximum rank estimation with two variables in Table 5 
         (39 hours for 200 point-estimates, 38 hours for CI bootstrapping 200)
    (3). matching maximum rank estimation with two variables for main 12 firms in Table 7
         (1 hours for 200 point-estimates)
  All results are restored in julia_merger_result folder.
  If you want to run the each estimation block, please switch
    want_to_run_plotting = "not_run" to "run"
    want_to_run_multivariate_estimation = "not_run" to "run"
    want_to_run_two_variable_estimation = "not_run" to "run"
    want_to_run_only_main_firms = "not_run" to "run"
  for plotting the shape of the objective function and (1), (2), (3).
  If you do not switch (= default), ship_merger_score_estimation.jl just reads current results 
  from julia_merger_results folder and generates (and updates) latex tables in ship_merger_table folder.
  
  This also generates Table 5, Table 7, and Table 10 in Appendix as Latex format in julia_merger_table folder.
  The script also translates results of Table 5 into comparable results as Latex format for Table 6.
  If you switch 
  temp_subsidy_type = "shared" into temp_subsidy_type = "to_buyer" on line 41,
  you can run all processes for robustness checks with another subsidy specification shown in Monte Carlo Section.

#-----------------------------------------#
# Counterfactual simulation (Julia)
#-----------------------------------------#
- shipmerger/ship_merger_counterfactual.jl
  It takes 100 hours (360402.808865 seconds) to compute 30 simulated matching outcome 
  for 54 scenarios (6 subsidy amounts * 9 subsidy thresholds). 
  
#--------------------------------------#
# Auxiliary files and datasets (Julia)
#--------------------------------------#
The following scripts are loaded only as module and data in ship_merger_score_estimation.jl.
- shipmerger/functions.jl
- shipmerger/functions_for_monte_carlo.jl (differs only in cost specification)
- shipmerger/output/data_for_maximum_rank_estimation.csv

The following intermediate scripts are used only in ship_merger_score_estimation.jl.
- shipmerger/ship_merger_score_estimation_plotting.jl
- shipmerger/ship_merger_score_estimation_point_estimates_and_confidence_intervals.jl
- shipmerger/ship_merger_score_estimation_generate_latex.jl
- shipmerger/ship_merger_score_estimation_generate_latex_all_firms.jl
- shipmerger/ship_merger_score_estimation_generate_latex_only_main_firms.jl


#-------------------------------#
# For Appendices (Julia)
#-------------------------------#
You can run the following files for generating results in Appendices.

- shipmerger/ship_merger_monte_carlo.jl

- shipmerger/ship_merger_comparative_statics.jl

- shipmerger/ship_merger_comparative_statics_large_firms.jl








