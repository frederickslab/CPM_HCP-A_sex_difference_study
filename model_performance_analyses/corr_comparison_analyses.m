% run t-tests to compare median R values between sexes

%% setup variables
% adjust the variables below as needed!
n_runs = 100;
param_list = {'ravlt_L','ravlt_IR','facename'};
scan_type_list = {'rfMRI_REST1_AP','rfMRI_REST1_PA','rfMRI_REST2_AP','rfMRI_REST2_PA','tfMRI_CARIT','tfMRI_FACENAME','tfMRI_VISMOTOR'};

%% load in R values
load('../BIG_data_from_CPM_HCP-Aging/all_medianR.mat', 'all_R_arrays');

%% run t-tests and collect necessary values

for param = 1:length(param_list)
    load(sprintf('../BIG_data_from_CPM_HCP-Aging/%s_cpm_output.mat',char(param_list{param})),'cpm_output')
    ttest_output = struct();
    
    mean_F_values = [];
    mean_M_values = [];
    h_values = [];
    p_values = [];
    tstat_values = [];

    for i = 1:length(scan_type_list)
        % extract female and male subj median R arrays
        arr_F = cpm_output.F_cpm_output.corr_struct.(char(scan_type_list{i}))(1,:);
        arr_M = cpm_output.M_cpm_output.corr_struct.(char(scan_type_list{i}))(1,:);
        
        mean_F = mean(arr_F);
        mean_M = mean(arr_M);
        mean_F_values = [mean_F_values; mean_F];
        mean_M_values = [mean_M_values; mean_M];

        [h,p,ci,stats] = ttest2(arr_F,arr_M);
        h_values = [h_values; h];
        p_values = [p_values; p];
        tstat_values = [tstat_values; stats.tstat];
    end

    ttest_output = struct('h', h_values, 'p', p_values, 'mean_F', mean_F_values, 'mean_M', mean_M_values, 'tstat', tstat_values);
    all_R_arrays.(char(sprintf('%s_R',char(param_list(param))))).(char('FM_ttest_results')) = ttest_output;
end
