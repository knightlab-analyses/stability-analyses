# stability-analyses
This set of codes are used for reproducing all the simulation studies and experimental microbiome data applications in Stability manuscript.

I. General code:
	
code_method folder: contain codes to reproduce simulation results for continuous outcomes

	getStability.R: function to calculate Stability Index

	cv_method.R: code for 4 selected feature selection methods with user-defined parameter grids and cross-validations for parameter tuning when applied to continuous outcomes

	cv_method_binary_update.R: code for 4 selected feature selection methods with user-defined parameter grids and cross-validations for parameter tuning when applied to binary outcomes

	stab_data_applications.R: function to perform hypothesis testing using bootstrap for continuous outcomes

	stab_data_applications_binary.R: function to perform hypothesis testing using bootstrap for binary outcomes

	bootstrap_test_compLasso_rf.R: general functions for comparing feature selection methods using hypothesis testing based on bootstrap when applied to continuous outcomes

	bootstrap_test_compLasso_rf_binary.R: general functions for comparing feature selection methods using hypothesis testing based on bootstrap when applied to binary outcomes

	source code for compositional lasso (continuous outcome) is available at: https://www.math.pku.edu.cn/teachers/linw/software.html
	source code for compositional lasso (binary outcome) is available at: https://github.com/UVic-omics/Microbiome-Variable-Selection


II. Simulation part (within simulations folder):

sim_data_generation folder: contain codes to generate simulated data

	sim_dat_ind_toeplitz: code to generate simulated data with Independent and Toeplitz correlation designs
	sim_dat_block.R: code to generate simulated data with Block correlation design
	run_sim_data.sh: bash commands for running simulation data generation code on HPC

code_sim_cts folder: contain codes to reproduce simulation results for continuous outcomes 

	cv_sim_apply.R: general functions for applying selected feature selection methods to simulated data when applied to continuous outcomes

	1. compute Stability and MSE for different simulation scenarios
	ind_results.R: code for comparing 3 methods (lasso, elastic net, random forests) in simulated data with Independent design and continuous outcomes
	toe_results.R: code for comparing 3 methods (lasso, elastic net, random forests) in simulated data with Toeplitz design and continuous outcomes
	block_results.R: code for comparing 3 methods (lasso, elastic net, random forests) in simulated data with Block design and continuous outcomes
	CL_sim_apply.R: code for obtaining results for compositional lasso in all simulation correlation designs with continuous outcomes

	2. hypothesis testing with bootstrap for selected simulation scenarios
	boot_CL_testing.R: code for calculating bootstrapped confidence interval for compositional lasso method in simulated data with continous outcomes
	boot_RF_testing.R: code for calculating bootstrapped confidence interval for random forests method in simulated data with continous outcomes

	3. bash commands
	run_sim_cts.sh: bash commands for running simulation code for continous outcomes on HPC


code_sim_bin folder: contain codes to reproduce simulation results for binary outcomes

	cv_sim_apply_binary_update.R: general functions for applying selected feature selection methods to simulated data when applied to binary outcomes

	1. compute Stability and AUC for different simulation scenarios
	ind_results_binary_update.R: code for comparing all 4 methods in simulated data with Independent design and binary outcomes
	toe_results_binary_update.R: code for comparing all 4 methods in simulated data with Toeplitz design and binary outcomes
	block_results_binary_update.R: code for comparing all 4 methods in simulated data with Block design and binary outcomes

	2. hypothesis testing with bootstrap for selected simulation scenarios
	boot_sim_binary.R: code for calculating bootstrapped confidence interval for compositional lasso and random forests methods in simulated data with binary outcomes

	3. bash commands
	run_sim_bin.sh: bash commands for running simulation code for binary outcomes on HPC

notebooks_sim_cts folder: contain notebooks (R) to summarize simulation results for continuous outcome

notebooks_sim_bin folder: contain notebooks (R) to summarize simulation results for binary outcome	

results_summary_cts folder: contain outputs of tables from notebooks in notebooks_sim_cts folder

results_summary_bin folder: contain outputs of tables from notebooks in notebooks_sim_bin folder

figures_combined folder: contain figures generated for both continous and binary outcomes based on notebook 6_make_figures_combined in notebooks_sim_bin folder


III. application part (within data_application folder):

	code_cts folder: contain code for real data applications to BMI & soil datasets for continuous outcomes

	code_bin folder: contain code for real data applications to BMI & soil datasets for binary outcomes

	notebooks_applications folder: contain notebooks (R) to summarize microbiome application results for continuous and binary outcomes

	88soils folder: contain data and application results for soil datast

	BMI folder: contain data and application results for BMI datast

