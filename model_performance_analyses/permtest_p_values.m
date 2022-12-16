% find permutation testing p-values

% from corey: calculate the number of times your performance metric 
%   of interest in the permuted data (below, 'rho_null') is greater than 
%   or equal to your test statistic (in my case, I chose 'median_rho', 
%   the median of my CPM performance in my original (unshuffled) data.

% tmp = length(find(rho_null >=median_rho));
% pval = tmp/1000;


clear all

%% setup variables
% adjust the variables below as needed!
n_runs = 1000;
param_list = {'facename','ravlt_L','ravlt_IR'};
scan_type_list = {'rfMRI_REST1_AP','rfMRI_REST1_PA','rfMRI_REST2_AP','rfMRI_REST2_PA','tfMRI_CARIT','tfMRI_FACENAME','tfMRI_VISMOTOR'};
subj_group_list = {'all', 'F', 'M'};
covariates = {'no_cov', 'motion', 'age','motion_age','motion_age_interact'};
covariate_of_choice = covariates{5};

for i = 1:length(param_list)
    tic;
    % load cpm permutation testing results
    load(sprintf('../BIG_data_from_CPM_HCP-Aging/%s_cpm_permtest_output_1000run_%s.mat',char(param_list{i}),covariate_of_choice),'cpm_permtest_output')
    
    % load spearman's rho of every median-performing CPM model
    load('../BIG_data_from_CPM_HCP-Aging/unpermuted_corr_stats.mat');
    
    for j = 1:length(subj_group_list)
        
        perm_data_mat = zeros(7,1000);
        med_rho_arr = zeros(7,1);
        rho_count_arr = zeros(7,1);
        pval_arr = zeros(7,1);

        for k = 1:length(scan_type_list)
            perm_data = cpm_permtest_output.(char(sprintf('%s_cpm_permtest_output',subj_group_list{j}))).corr_struct.(char(scan_type_list{k}))(1,:); % permuted data (n_runs # of values)
            med_rho = corr_stats.(char(sprintf('%s_stats',param_list{i}))).(char(subj_group_list{j})).spearmans_median_coeff(k); % median of unpermuted data (single value)
            rho_count = length(find(perm_data >= med_rho)); % # of values in permuted data that are larger than median of unpermuted data
            pval = rho_count/n_runs;

            perm_data_mat(k,:) = perm_data;
            med_rho_arr(k) = med_rho;
            rho_count_arr(k) = rho_count;
            pval_arr(k) = pval;
        end
        permtest_p_struct.(char(sprintf('%s_permtest_p',param_list{i}))).(char(sprintf('%s_group',subj_group_list{j}))).(char('perm_data')) = perm_data_mat;
        permtest_p_struct.(char(sprintf('%s_permtest_p',param_list{i}))).(char(sprintf('%s_group',subj_group_list{j}))).(char('median_rho')) = med_rho_arr;
        permtest_p_struct.(char(sprintf('%s_permtest_p',param_list{i}))).(char(sprintf('%s_group',subj_group_list{j}))).(char('rho_count')) = rho_count_arr;
        permtest_p_struct.(char(sprintf('%s_permtest_p',param_list{i}))).(char(sprintf('%s_group',subj_group_list{j}))).(char('pval')) = pval_arr;
    end
    fprintf('Retrieved p value for %s\n',param_list{i})
    toc;
end



