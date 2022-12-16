close all
% clear
% clc

%% load the 268 atlas based connectivity matrices
data_path = '../../BIG_data_from_CPM_HCP-Aging/all_conn_mats/';
load([data_path, 'all_conn_mats_indiv.mat'], 'conn_mat_struct_all');
load([data_path, 'pt_id_connmat_grouped.mat'], 'pt_id_connmat_struct');

scan_type_list = {'rfMRI_REST1_AP', 'rfMRI_REST1_PA', 'rfMRI_REST2_AP', 'rfMRI_REST2_PA','tfMRI_CARIT', 'tfMRI_FACENAME', 'tfMRI_VISMOTOR'};

for s = 1:length(scan_type_list)
    [mat_all, mat_F, mat_M] = concat_connmats(conn_mat_struct_all, pt_id_connmat_struct, scan_type_list{s});
    
    % compute association ratios for each group
    a_ratio_all = compute_association_ratio(mat_all);
    a_ratio_F = compute_association_ratio(mat_F);
    a_ratio_M = compute_association_ratio(mat_M);
    
%     % plot association ratios of each group
%     plot_association_ratio_hm(a_ratio_all, scan_type_list{s}, s, 'all');
%     plot_association_ratio_hm(a_ratio_F, scan_type_list{s}, s,'F');
%     plot_association_ratio_hm(a_ratio_M, scan_type_list{s}, s,'M');
    
    % calculate mean of association ratios for each network in sex-based groups
    a_ratio_F_mean = mean(a_ratio_F);
    a_ratio_M_mean = mean(a_ratio_M);

    % perform ttest to compare association ratios of networks between sexes for each scan type
    [h, p, ci, stats] = ttest2(a_ratio_F, a_ratio_M);
    ttest_output = struct('h', h, 'p', p, 'ci', ci, 'mean_F', a_ratio_F_mean, 'mean_M', a_ratio_M_mean, 'stats', stats);
    a_ratio_ttest_results.(char(scan_type_list(s))) = ttest_output; % outputs are all 1x10 vectors, where each element is a separate network
    
    for i=1:10
        [p, h, stats] = ranksum(a_ratio_F(:,i), a_ratio_M(:,i));
        rs_p(i) = p;
        rs_h(i) = h;
        rs_stats(i) = stats;
        rs_output = struct('p', rs_p, 'h', rs_h, 'stats', rs_stats);
        a_ratio_rs_results.(char(scan_type_list(s))) = rs_output;
    end
end

%% fxn to concatenate all existing conn mats for each input scan type name
function [all_mat, F_mat, M_mat] = concat_connmats(connmat_struct, pt_id_struct, scan_type_name)
    % get all pt id's for each subj group
    fn_all = pt_id_struct.(char(scan_type_name)).all_subs;
    fn_F = pt_id_struct.(char(scan_type_name)).F_subs;
    fn_M = pt_id_struct.(char(scan_type_name)).M_subs;
    
    % initialize matrices for each group
    all_mat = zeros(268,268,size(fn_all,1));
    F_mat = zeros(268,268,size(fn_F,1));
    M_mat = zeros(268,268,size(fn_M,1));
    
    % fill in matrices for each group
    for s = 1:length(fn_all)
        all_mat(:,:,s) = connmat_struct.(char(scan_type_name)).(char(fn_all{s}));
    end
    for s = 1:length(fn_F)
        F_mat(:,:,s) = connmat_struct.(char(scan_type_name)).(char(fn_F{s}));
    end
    for s = 1:length(fn_M)
        M_mat(:,:,s) = connmat_struct.(char(scan_type_name)).(char(fn_M{s}));
    end
end

%% fxn to plot heatmap of association ratios across all networks for each scan type
function plot_association_ratio_hm(A_ratio, scan_type_name, scan_type_num, group_name)
    switch scan_type_num
            case 1
                scan_type = 'rfMRI\_REST1\_AP';
            case 2
                scan_type = 'rfMRI\_REST1\_PA';
            case 3
                scan_type = 'rfMRI\_REST2\_AP';
            case 4
                scan_type = 'rfMRI\_REST2\_PA';
            case 5
                scan_type = 'tfMRI\_CARIT';
            case 6
                scan_type = 'tfMRI\_FACENAME';
            case 7
                scan_type = 'tfMRI\_VISMOTOR';
    end

    fig = figure;
    set(gcf, 'Position', [0, 0, 500, 1500]);
    
    hm = heatmap(A_ratio,'Colormap', jet);
    hm.Title = [scan_type];
    hm.XDisplayLabels = {'MF','FP','DMN','Mot','VI','VII','VAs','SAL','SC','Cer'};
    hm.XLabel = '10 networks';
    hm.YLabel = 'subject number';
    
    hm_filename = sprintf('../cpm_figures/network_segregation_heatmaps/netseg_hm_%s_%s.png', scan_type_name, group_name);
    saveas(fig,hm_filename) 
end

%% fxn to compute association ratio
function ass_ratio = compute_association_ratio(mat_all)
    % set up variables for number of subjects, nodes, and networks
    no_sub = size( mat_all, 3);
    no_nodes = size( mat_all, 1);
    no_net = 10;
    
    % load the network labels for the 268 parcellation
    % first column is the ROI index
    % second column is the 10 network labels (full network names commented below)
    % convert the two column lookup table to 1D label
    %     1 = medial frontal
    %     2 = frontoparietal
    %     3 = default mode
    %     4 = motor cortex, 
    %     5 = visual A
    %     6 = visual B
    %     7 = visual association
    %     8 = salience
    %     9 = subcortical
    %     10 = cerebellum
    look_table = dlmread('Shen268_10network');
    
    % set sigma ***SUYEON, FIGURE OUT WHAT THIS IS EXACTLY...
    sigma = mean( 1- mat_all(1:50:end));
    
    % 
    label = zeros(1, no_nodes);
    for i = 1:no_net
        cur_id= find( look_table(:,2) ==i);
        label( look_table(cur_id,1)) = i;
    end
    
    % initialize association ratio matrix (no_sub x no_net)
    ass_ratio = zeros(no_sub, no_net);
    
    % fill in association ratio matrix
    for i = 1:no_sub
        cur_mat = squeeze( mat_all(:,:,i));
        ass_ratio(i, :) = association( cur_mat, sigma, label);
    end
end
