# public
# Estimating Endogenous Coalitional Mergers: Merger Costs and Assortativeness of Size and Specialization.
# written by Suguru Otani, 2021.08.05

PATH must be on shipmerger/.
The following folders contain all R implementations and raw data.
- shipmerger/main: This contains execution files.
- shipmerger/input: This contains a raw dataset.
- shipmerger/output: This contains cleaned dataset for Julia.
- shipmerger/figuretable: This contains figures and tables shown in my paper.
- shipmerger/figure_table_ppt: This contains power point file creating figures of sankey diagrams.

Since computing all processes of Julia on a single laptop may take a few weeks,
I locate the following folders restoring current Julia result files on ../shipmerger/.
- shipmerger/julia_merger_result : This contains all intermediate files.
- shipmerger/julia_merger_table: This contains all tables converting intermediate files into latex style.
                                 XXX_generate_latex_YYY named functions are only used for the conversion.
- shipmerger/julia_merger_figure: This contains all figures.

The following files restore (and update current files if you run a code) the results on these folders. 
You can replicate all results by following processes sequentially.
However, as randomness of DE algorithm and sampling, some results may vary. 
To keep a record of my final results, I restore all intermediate files in shipmerger/julia_merger_result.


## Data Cleaning and Reduced Form Analysis (Rproject)
- shipmerger/main/01_make_figuretable.R
  This cleans raw data file (210627firmdata_processed.csv) in shipmerger/input.
  Then, this produces cleaned data file (data_for_maximum_rank_estimation.csv) for Julia in shipmerger/output.
  This also generates related Figures and Tables in my paper on shipmerger/figuretable.
  It takes just a minute.
- shipmerger/main/02_make_sankey_diagram
  This generates sankey diagram in Figure 6. 
  Input data is manually assigned based on results on ship_merger_counterfactual.jl.


## Matching Maximum Rank Estimation (Julia)

- ship_merger_score_estimation_HHI.jl
    (1). the simplest matching maximum rank estimation in Appendix
         (4 hours for 50 point-estimates, 6 hours for CI bootstrapping 200)
- shipmerger/ship_merger_score_estimation.jl
  This computes 8 specifications with bootstrap 95% CI for
    (1). multivariate matching maximum rank estimation in Appendix
         (22 hours for 200 point-estimates, 21 hours for CI bootstrapping 200)
    (2). matching maximum rank estimation with two variables in Table 5 
         (39 hours for 200 point-estimates, 38 hours for CI bootstrapping 200)
    (3). matching maximum rank estimation with two variables for main 12 firms in Table 7
         (1 hours for 200 point-estimates)
  All results are restored in julia_merger_result folder.
  If you want to run the each estimation block, please switch
    want_to_run_multivariate_estimation = "not_run" to "run"
    want_to_run_two_variable_estimation = "not_run" to "run"
    want_to_run_only_main_firms = "not_run" to "run"
  If you do not switch (= default), ship_merger_score_estimation.jl just reads current results 
  from julia_merger_results folder and generates (and updates) latex tables in ship_merger_table folder.
  This generates Table 4, 5, 6 in the main text and Table 11, 12 in Appendix 
  as Latex format in julia_merger_table folder.
  If you switch 
  temp_subsidy_type = "shared" into temp_subsidy_type = "to_buyer" on line 40,
  you can run all processes for robustness checks with another subsidy specification shown in Monte Carlo Section.


## Counterfactual simulation (Julia)
- shipmerger/ship_merger_counterfactual.jl
  It takes 100 hours (360402.808865 seconds) to compute 20 simulated matching outcome 
  for 60 scenarios (6 subsidy amounts * 10 subsidy thresholds). 




## Auxiliary files and datasets (Julia)

The following scripts are loaded only as module and data in ship_merger_score_estimation.jl.
- shipmerger/functions.jl
- shipmerger/output/data_for_maximum_rank_estimation.csv

The following intermediate scripts are used only in ship_merger_score_estimation.jl.
- shipmerger/ship_merger_score_estimation_plotting.jl
- shipmerger/ship_merger_score_estimation_point_estimates_and_confidence_intervals.jl
- shipmerger/ship_merger_score_estimation_generate_latex.jl
- shipmerger/ship_merger_score_estimation_generate_latex_all_firms.jl
- shipmerger/ship_merger_score_estimation_generate_latex_only_main_firms.jl

The following intermediate scripts are used only in shipmerger/ship_merger_counterfactual.jl.
- shipmerger/ship_merger_counterfactual_plotting.jl
- shipmerger/ship_merger_counterfactual_for_sankey_diagrams.jl
- shipmerger/ship_merger_counterfactual_functions.jl


# For Appendices (Julia)
#-------------------------------#
You can run the following files for generating results in Appendices.
- shipmerger/ship_merger_score_estimation_HHI.jl
- shipmerger/ship_merger_monte_carlo.jl
- shipmerger/ship_merger_comparative_statics.jl

- shipmerger/ship_merger_monte_carlo_large_firms
- shipmerger/ship_merger_comparative_statics_large_firms.jl

The following intermediate scripts are used only in the above files.
- shipmerger/ship_merger_score_estimation_point_estimates_and_confidence_intervals_HHI
- shipmerger/functions_for_monte_carlo.jl (differs only in cost specification)
- shipmerger/ship_merger_counterfactual_functions.jl
