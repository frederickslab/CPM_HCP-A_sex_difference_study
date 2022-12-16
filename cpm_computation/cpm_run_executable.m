% executable file for all scripts to run in MRRC server

clear all

% all options:
run_all_conn_mats_bool = 0; % includes 'get_all_conn_mats' fxn
cpm_run_bool = 1;
cpm_permutation_testing_bool = 0;

% all arguments:
scan_type_list = {'rfMRI_REST1_AP','rfMRI_REST1_PA','rfMRI_REST2_AP','rfMRI_REST2_PA','tfMRI_CARIT','tfMRI_FACENAME','tfMRI_VISMOTOR'};
% param_list = {'ravlt_L','ravlt_IR','facename','neon'};
param_list = {'facename'};
var_list = {'all', 'sex'};
hcp_a_or_a4 = 'hcp_a';
n_runs = 1000;
p_thresh = 0.01;
k_folds = 5;
corr_type = 'partial';

% all functions:
if run_all_conn_mats_bool
    run_all_conn_mats(scan_type_list, hcp_a_or_a4);
end

if cpm_run_bool
    cpm_output = cpm_run(param_list, scan_type_list, var_list, hcp_a_or_a4, n_runs, p_thresh, k_folds, corr_type);
end

if cpm_permutation_testing_bool
    cpm_permtest_output = cpm_permutation_testing(param_list, scan_type_list, var_list, hcp_a_or_a4, n_runs, p_thresh, k_folds, corr_type);
end
