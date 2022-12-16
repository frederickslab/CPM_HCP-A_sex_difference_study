% visualize pmasks of significant edges from cpm

% outputs:
%   saves 10-network consensus heatmaps of positive and negative matrices for each scan type
%   of each parameter in `param_list` (for later use in BIS Connviewer)

%% Implementation

close all
clear all
% clearvars -except cpm_output

param_list = {'facename','ravlt_L','ravlt_IR'};
scan_type_list = {'rfMRI_REST1_AP', 'rfMRI_REST1_PA', 'rfMRI_REST2_AP', 'rfMRI_REST2_PA','tfMRI_CARIT', 'tfMRI_FACENAME', 'tfMRI_VISMOTOR'};
covariates = {'no_cov', 'motion', 'age','motion_age','motion_age_interact'};
covariate_of_choice = covariates{5};

p_thresh = 0.01; % p threshold used for CPM run
k_folds = 5; % number of k folds used in CPM run
n_runs = 1000; % number of CPM runs
k_fold_threshold = 2; % threshold of number of edge selections out of k_folds folds (ie, 3 selections out of 5 folds)
iter_threshold = 400; % threshold of number of edge selections out of trial_count CPM iterations/runs (ie, 600 out of 1000 runs)

standardized_or_nah = 0; % 0 = colorbars NOT standardized between all figs; 1 = colorbars ARE standardized between all figs
plot_pmask_separately_bool = 0;
plot_F_M_difference_pmasks_bool = 1;

for n = 1:length(param_list)
    tic;

    load(sprintf('../BIG_data_from_CPM_HCP-Aging/pmask_structs/%s_pmask_1000run_%s.mat',param_list{n},covariate_of_choice))
    
    %% set all necessary variables
    param = char(param_list{n}); % char variable with name of param
        
    %% get positive and negative pmasks and their sizes
    for i = 6 %1:length(scan_type_list)
        % get pos and neg matrices from F and M models
        [pos_mat_F,neg_mat_F,pos_mat_size_F,neg_mat_size_F] = get_consensus_mask(pmasks.F_pmasks.(char(scan_type_list{i})),k_folds,n_runs,k_fold_threshold,iter_threshold);
        [pos_mat_M,neg_mat_M,pos_mat_size_M,neg_mat_size_M] = get_consensus_mask(pmasks.M_pmasks.(char(scan_type_list{i})),k_folds,n_runs,k_fold_threshold,iter_threshold);
        
        % print number of sig edges in pos and neg matrices
        pos_mat_size_F
        neg_mat_size_F      
        pos_mat_size_M
        neg_mat_size_M

        % plot pmask heatmaps separately for F and M
        if plot_pmask_separately_bool
            % plot F pmasks
            pmask_visualization(pos_mat_F,neg_mat_F, param, scan_type_list{i},'F', covariate_of_choice,standardized_or_nah);
            % plot M pmasks
            pmask_visualization(pos_mat_M,neg_mat_M, param, scan_type_list{i},'M', covariate_of_choice,standardized_or_nah);
        end
        
        % plot F-M difference pmasks
        if plot_F_M_difference_pmasks_bool
            % plot pos mat (F-M difference)
            pmask_difference_visualization(pos_mat_F,pos_mat_M,param,scan_type_list{i},covariate_of_choice,standardized_or_nah, 'positive');
            % plot neg mat (F-M difference)
            pmask_difference_visualization(neg_mat_F,neg_mat_M,param,scan_type_list{i},covariate_of_choice,standardized_or_nah, 'negative');
        end
    end
    fprintf('Retrieved pmask sizes for %s\n',param_list{n})
    toc;
end