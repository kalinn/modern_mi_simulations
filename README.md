# modern_mi_simulations
Analytic code used to produce results in the paper, "Performance of multiple imputation using modern machine learning methods in electronic health records data"

Here's an outline of the steps to reproduce the results:
1) Run the scripts in folder data_generation to create full plasmode simulation datasets and datasets with MCAR, MAR, and MNAR missingness.
2) Implement MICE, random forest multiple imputation, and DAE multiple imputation, followed by downstream Cox proportional hazards modeling and application of Rubin's rules using scripts in method_implementation. Running each script once will produce results from one simulation iteration. We ran each script 1,000 times in parallel on a high performance computing cluster.
3) Run results/ConsolidateResults.R to reduce information from simulation iterations into final data frames containing summary information.
4) Create plots and tables using results/MIPaperResults.Rmd
5) Table 1 was created using the script table1/table1.R
