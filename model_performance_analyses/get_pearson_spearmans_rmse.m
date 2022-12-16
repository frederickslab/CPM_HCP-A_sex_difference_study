% calculate pearson, spearman's, and rmse stats
clear all

%% setup variables
% adjust the variables below as needed!
n_runs = 1000;
param_list = {'facename','ravlt_L','ravlt_IR'};
scan_type_list = {'rfMRI_REST1_AP','rfMRI_REST1_PA','rfMRI_REST2_AP','rfMRI_REST2_PA','tfMRI_CARIT','tfMRI_FACENAME','tfMRI_VISMOTOR'};
covariates = {'no_cov', 'motion', 'age','motion_age','motion_age_interact'};
covariate_of_choice = covariates{5};

%% booleans to choose which fxns you'd like to run
permtest_or_nah = 1; % 0 = regular cpm, 1 = cpm permutation testing
derive_pearson_spearmans_rmse_stats_bool = 0;
plot_median_performances_bool = 1;
wilcoxon_rs_test_bool = 0;

%% collect pearson, spearman's, and rmse stats for all predictive models
if derive_pearson_spearmans_rmse_stats_bool
    for i = 1:length(param_list)
        tic;
    
        all_behav_score_table = readtable(sprintf('../HCP-A_info/memory_scores/%s_all.csv', param_list{i}));
        F_behav_score_table = readtable(sprintf('../HCP-A_info/memory_scores/%s_F.csv', param_list{i}));
        M_behav_score_table = readtable(sprintf('../HCP-A_info/memory_scores/%s_M.csv', param_list{i}));
        
    %     behav_tables = struct('all_behav_score_table',all_behav_score_table,'F_behav_score_table',F_behav_score_table,'M_behav_score_table',M_behav_score_table);
    
        % reg cpm
        if ~permtest_or_nah
            load(sprintf('../BIG_data_from_CPM_HCP-Aging/cpm_output/%s_cpm_output_1000run_%s.mat',char(param_list{i}),covariate_of_choice),'cpm_output')
    
            for j = 1:length(scan_type_list)
                %all
                actual_y = all_behav_score_table.(char(sprintf('%s_score',param_list{i})));
                predicted_y = cpm_output.all_cpm_output.y_hat_struct.(char(scan_type_list{j}));
                [all_corr_stats.(char(scan_type_list{j})).pearson_stats,all_corr_stats.(char(scan_type_list{j})).spearmans_stats,all_corr_stats.(char(scan_type_list{j})).rmse_stats] = pearson_spearman_rmse_calc(actual_y, predicted_y, n_runs);
                [all_corr_stats.pearson_median_coeff(j),all_corr_stats.pearson_median_p(j),all_corr_stats.pearson_median_idx(j),all_corr_stats.pearson_std_dev(j)] = get_median_stats(all_corr_stats.(char(scan_type_list{j})).pearson_stats, n_runs);
                [all_corr_stats.spearmans_median_coeff(j),all_corr_stats.spearmans_median_p(j),all_corr_stats.spearmans_median_idx(j),all_corr_stats.spearmans_std_dev(j)] = get_median_stats(all_corr_stats.(char(scan_type_list{j})).spearmans_stats, n_runs);
                rmse_sorted = sort(all_corr_stats.(char(scan_type_list{j})).rmse_stats);
                all_corr_stats.rmse_median(j) = rmse_sorted((n_runs/2) + 1);
    
                %F
                actual_y = F_behav_score_table.(char(sprintf('%s_score',param_list{i})));
                predicted_y = cpm_output.F_cpm_output.y_hat_struct.(char(scan_type_list{j}));
                [F_corr_stats.(char(scan_type_list{j})).pearson_stats,F_corr_stats.(char(scan_type_list{j})).spearmans_stats,F_corr_stats.(char(scan_type_list{j})).rmse_stats] = pearson_spearman_rmse_calc(actual_y, predicted_y, n_runs);
                [F_corr_stats.pearson_median_coeff(j),F_corr_stats.pearson_median_p(j),F_corr_stats.pearson_median_idx(j),F_corr_stats.pearson_std_dev(j)] = get_median_stats(F_corr_stats.(char(scan_type_list{j})).pearson_stats, n_runs);
                [F_corr_stats.spearmans_median_coeff(j),F_corr_stats.spearmans_median_p(j),F_corr_stats.spearmans_median_idx(j),F_corr_stats.spearmans_std_dev(j)] = get_median_stats(F_corr_stats.(char(scan_type_list{j})).spearmans_stats, n_runs);
                rmse_sorted = sort(F_corr_stats.(char(scan_type_list{j})).rmse_stats);
                F_corr_stats.rmse_median(j) = rmse_sorted((n_runs/2) + 1);
    
                %M
                actual_y = M_behav_score_table.(char(sprintf('%s_score',param_list{i})));
                predicted_y = cpm_output.M_cpm_output.y_hat_struct.(char(scan_type_list{j}));
                [M_corr_stats.(char(scan_type_list{j})).pearson_stats,M_corr_stats.(char(scan_type_list{j})).spearmans_stats,M_corr_stats.(char(scan_type_list{j})).rmse_stats] = pearson_spearman_rmse_calc(actual_y, predicted_y, n_runs);
                [M_corr_stats.pearson_median_coeff(j),M_corr_stats.pearson_median_p(j),M_corr_stats.pearson_median_idx(j),M_corr_stats.pearson_std_dev(j)] = get_median_stats(M_corr_stats.(char(scan_type_list{j})).pearson_stats, n_runs);
                [M_corr_stats.spearmans_median_coeff(j),M_corr_stats.spearmans_median_p(j),M_corr_stats.spearmans_median_idx(j),M_corr_stats.spearmans_std_dev(j)] = get_median_stats(M_corr_stats.(char(scan_type_list{j})).spearmans_stats, n_runs);
                rmse_sorted = sort(M_corr_stats.(char(scan_type_list{j})).rmse_stats);
                M_corr_stats.rmse_median(j) = rmse_sorted((n_runs/2) + 1);
    
            end
    
            corr_stats.(char(sprintf('%s_stats', param_list{i}))) = struct('all', all_corr_stats, 'F', F_corr_stats, 'M', M_corr_stats);
            fprintf('Collected all Pearson, Spearmans, and RMSE stats for %s\n', param_list{i})
        end
    
        % permtest
        if permtest_or_nah
            load(sprintf('../BIG_data_from_CPM_HCP-Aging/cpm_permtest_output/%s_cpm_permtest_output_1000run_%s.mat',char(param_list{i}),covariate_of_choice),'cpm_permtest_output')
    
            for j = 1:length(scan_type_list)
                %all
                actual_y = all_behav_score_table.(char(sprintf('%s_score',param_list{i})));
                predicted_y = cpm_permtest_output.all_cpm_permtest_output.y_hat_struct.(char(scan_type_list{j}));
                [all_corr_stats.(char(scan_type_list{j})).pearson_stats,all_corr_stats.(char(scan_type_list{j})).spearmans_stats,all_corr_stats.(char(scan_type_list{j})).rmse_stats] = pearson_spearman_rmse_calc(actual_y, predicted_y, n_runs);
                [all_corr_stats.pearson_median_coeff(j),all_corr_stats.pearson_median_p(j),all_corr_stats.pearson_median_idx(j),all_corr_stats.pearson_std_dev(j)] = get_median_stats(all_corr_stats.(char(scan_type_list{j})).pearson_stats, n_runs);
                [all_corr_stats.spearmans_median_coeff(j),all_corr_stats.spearmans_median_p(j),all_corr_stats.spearmans_median_idx(j),all_corr_stats.spearmans_std_dev(j)] = get_median_stats(all_corr_stats.(char(scan_type_list{j})).spearmans_stats, n_runs);
                rmse_sorted = sort(all_corr_stats.(char(scan_type_list{j})).rmse_stats);
                all_corr_stats.rmse_median(j) = rmse_sorted((n_runs/2) + 1);
    
                %F
                actual_y = F_behav_score_table.(char(sprintf('%s_score',param_list{i})));
                predicted_y = cpm_permtest_output.F_cpm_permtest_output.y_hat_struct.(char(scan_type_list{j}));
                [F_corr_stats.(char(scan_type_list{j})).pearson_stats,F_corr_stats.(char(scan_type_list{j})).spearmans_stats,F_corr_stats.(char(scan_type_list{j})).rmse_stats] = pearson_spearman_rmse_calc(actual_y, predicted_y, n_runs);
                [F_corr_stats.pearson_median_coeff(j),F_corr_stats.pearson_median_p(j),F_corr_stats.pearson_median_idx(j),F_corr_stats.pearson_std_dev(j)] = get_median_stats(F_corr_stats.(char(scan_type_list{j})).pearson_stats, n_runs);
                [F_corr_stats.spearmans_median_coeff(j),F_corr_stats.spearmans_median_p(j),F_corr_stats.spearmans_median_idx(j),F_corr_stats.spearmans_std_dev(j)] = get_median_stats(F_corr_stats.(char(scan_type_list{j})).spearmans_stats, n_runs);
                rmse_sorted = sort(F_corr_stats.(char(scan_type_list{j})).rmse_stats);
                F_corr_stats.rmse_median(j) = rmse_sorted((n_runs/2) + 1);
    
                %M
                actual_y = M_behav_score_table.(char(sprintf('%s_score',param_list{i})));
                predicted_y = cpm_permtest_output.M_cpm_permtest_output.y_hat_struct.(char(scan_type_list{j}));
                [M_corr_stats.(char(scan_type_list{j})).pearson_stats,M_corr_stats.(char(scan_type_list{j})).spearmans_stats,M_corr_stats.(char(scan_type_list{j})).rmse_stats] = pearson_spearman_rmse_calc(actual_y, predicted_y, n_runs);
                [M_corr_stats.pearson_median_coeff(j),M_corr_stats.pearson_median_p(j),M_corr_stats.pearson_median_idx(j),M_corr_stats.pearson_std_dev(j)] = get_median_stats(M_corr_stats.(char(scan_type_list{j})).pearson_stats, n_runs);
                [M_corr_stats.spearmans_median_coeff(j),M_corr_stats.spearmans_median_p(j),M_corr_stats.spearmans_median_idx(j),M_corr_stats.spearmans_std_dev(j)] = get_median_stats(M_corr_stats.(char(scan_type_list{j})).spearmans_stats, n_runs);
                rmse_sorted = sort(M_corr_stats.(char(scan_type_list{j})).rmse_stats);
                M_corr_stats.rmse_median(j) = rmse_sorted((n_runs/2) + 1);
            end
    
            corr_stats.(char(sprintf('%s_stats', param_list{i}))) = struct('all', all_corr_stats, 'F', F_corr_stats, 'M', M_corr_stats);
            fprintf('Collected all Pearson, Spearmans, and RMSE stats (from permutation testing) for %s\n', param_list{i})
        end
    
        toc;
    end
    
    if ~permtest_or_nah
        save('../cpm_analyses/model_performance_stats_analyses/unpermuted_corr_stats.mat','corr_stats');
        fprintf('Saved unpermuted data correlation stats!\n')
    else
        save('../cpm_analyses/model_performance_stats_analyses/permuted_corr_stats.mat','corr_stats');
        fprintf('Saved permuted data correlation stats!\n')
    end
end

%% plot Pearson, Spearman's, and RMSE stats
if plot_median_performances_bool
    close all

    if ~permtest_or_nah
        load('../cpm_analyses/model_performance_stats_analyses/unpermuted_corr_stats.mat');
    else
        load('../cpm_analyses/model_performance_stats_analyses/permuted_corr_stats.mat');
    end
    plot_median_performances(corr_stats, 'pearson', scan_type_list, permtest_or_nah)
    plot_median_performances(corr_stats, 'spearmans', scan_type_list, permtest_or_nah)
end

%% collect pearson, spearman's, and rmse stats for all predictive models
if wilcoxon_rs_test_bool
    load('../cpm_analyses/model_performance_stats_analyses/unpermuted_corr_stats.mat');
    for i = 1:length(param_list)
        tic;
        for j = 1:length(scan_type_list)
            [wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_p')(j),wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_z')(j),wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_rs')(j)] = wilcoxon_rs_test(corr_stats.(char(sprintf('%s_stats',param_list{i}))).F.(char(scan_type_list{j})).pearson_stats, corr_stats.(char(sprintf('%s_stats',param_list{i}))).M.(char(scan_type_list{j})).pearson_stats, n_runs);
            [wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_p')(j),wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_z')(j),wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_rs')(j)] = wilcoxon_rs_test(corr_stats.(char(sprintf('%s_stats',param_list{i}))).F.(char(scan_type_list{j})).spearmans_stats, corr_stats.(char(sprintf('%s_stats',param_list{i}))).M.(char(scan_type_list{j})).spearmans_stats, n_runs);
        end
        wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_p') = wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_p')';
        wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_z') = wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_z')';
        wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_rs') = wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('pearson_rs')';
        wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_p') = wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_p')';
        wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_z') = wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_z')';
        wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_rs') = wilcoxon_p_values.(char(sprintf('%s',param_list{i}))).('spearmans_rs')';
        fprintf('Ran nonparametric t-test comparing F and M correlations for %s\n', param_list{i})
        toc;
    end
    save('../cpm_analyses/model_performance_stats_analyses/wilcoxon_FM_pvalues.mat', 'wilcoxon_p_values')
end
%% fxn to calculate all Pearson, Spearman's, and RMSE stats (comparing y-observed to y-hat)
function [pearson_stats,spearmans_stats,rmse_stats] = pearson_spearman_rmse_calc(y_observed, y_hat, n_iter)
%y_observed = M x 1 matrix, where M is n_subjs
%y_hat = M x N matrix, where M is n_subjs and N is number of CPM iterations (in our case, 1000)

pearson_stats = zeros(2, n_iter);
spearmans_stats = zeros(2, n_iter);
rmse_stats = zeros(1, n_iter);

for i = 1:n_iter
    %calc pearson R and p
    [pearson_stats(1,i),pearson_stats(2,i)] = corr(y_hat(:,i),y_observed);
    
    %calc spearman's rho and p
    [spearmans_stats(1,i),spearmans_stats(2,i)] = corr(y_hat(:,i),y_observed,'Type','Spearman');
    
    %calc rmse
    rmse_stats(1,i) = sqrt(mean((y_observed - y_hat(:,i)).^2));
end

end

%% fxn to get Pearson and Spearman's stats only for median-performing models
function [median_coeff, median_p, median_coeff_idx, std_dev] = get_median_stats(stats_arr,n_runs)

corr_sorted = sort(stats_arr(1,:));
median_coeff = corr_sorted((n_runs/2) + 1);
median_coeff_idx = find(stats_arr(1,:) == median_coeff);
median_p = stats_arr(2,median_coeff_idx);
std_dev = std(corr_sorted);

end

%% fxn to plot median-performing model stats data (across all scan types)
function plot_median_performances(corr_stats_struct, stats_type, scan_type_list, permtest_or_nah)
fig_position_size = [1000, 500, 1200, 400];
fs = 12; %fontsize
title_fs = 12;

if strcmp(stats_type, 'pearson')
    %% whole group 
    ylims = [-0.05 0.45];

    fig = figure();
    set(gcf, 'Position', fig_position_size,'color','w'); % fullforpptsize
    t = tiledlayout(1,3,'TileSpacing','compact');
    
    % facename subplot
    nexttile
    data = corr_stats_struct.facename_stats.all.pearson_median_coeff';
    errlow = corr_stats_struct.facename_stats.all.pearson_std_dev';
    errhigh = corr_stats_struct.facename_stats.all.pearson_std_dev';
    
    b4 = bar(data, 'facecolor', [0.4940 0.1840 0.5560]);
    ax4 = gca;
    set(ax4,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    if ~permtest_or_nah
        ylim(ylims);
    end
    title('Pearson R of Median-Performing FN-TR Model','FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    hold on
    er = errorbar(1:numel(scan_type_list),data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    
    % ravlt_L subplot
    nexttile
    data = corr_stats_struct.ravlt_L_stats.all.pearson_median_coeff';
    errlow = corr_stats_struct.ravlt_L_stats.all.pearson_std_dev';
    errhigh = corr_stats_struct.ravlt_L_stats.all.pearson_std_dev';

    b1 = bar(data, 'facecolor', [0.4940 0.1840 0.5560]);
    ax1 = gca;
    set(ax1,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title('Pearson R of Median-Performing RAVLT-L Model','FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    hold on
    er = errorbar(1:numel(scan_type_list),data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    
    % ravlt_IR subplot
    nexttile
    data = corr_stats_struct.ravlt_IR_stats.all.pearson_median_coeff';
    errlow = corr_stats_struct.ravlt_IR_stats.all.pearson_std_dev';
    errhigh = corr_stats_struct.ravlt_IR_stats.all.pearson_std_dev';

    b2 = bar(data, 'facecolor', [0.4940 0.1840 0.5560]);
    ax2 = gca;
    set(ax2,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title('Pearson R of Median-Performing RAVLT-IR Model','FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    hold on
    er = errorbar(1:numel(scan_type_list),data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    
    % link y-axes of all subplots
    linkaxes([ax4 ax1 ax2], 'y');
    
    % save figure
    if ~permtest_or_nah
        saveas(fig,'../cpm_figures/R_comparison_figs/median/median_pearson_comparison_all_subjs_wErrorbars.png')
    end
    
    if permtest_or_nah
        saveas(fig,'../cpm_figures/R_comparison_figs/median/median_pearson_comparison_all_subjs_permtest_wErrorbars.png')
    end
    
    %% sex-based group
    ylims = [-0.1 0.5];
    
    % combine female and male arrays to create stacked bar graphs
    ravlt_L_FM_arr = cat(2, corr_stats_struct.ravlt_L_stats.F.pearson_median_coeff',corr_stats_struct.ravlt_L_stats.M.pearson_median_coeff');
    ravlt_IR_FM_arr = cat(2, corr_stats_struct.ravlt_IR_stats.F.pearson_median_coeff',corr_stats_struct.ravlt_IR_stats.M.pearson_median_coeff');
    facename_FM_arr = cat(2, corr_stats_struct.facename_stats.F.pearson_median_coeff',corr_stats_struct.facename_stats.M.pearson_median_coeff');
    
    ravlt_L_FM_err_arr = cat(2, corr_stats_struct.ravlt_L_stats.F.pearson_std_dev',corr_stats_struct.ravlt_L_stats.M.pearson_std_dev');
    ravlt_IR_FM_err_arr = cat(2, corr_stats_struct.ravlt_IR_stats.F.pearson_std_dev',corr_stats_struct.ravlt_IR_stats.M.pearson_std_dev');
    facename_FM_err_arr = cat(2, corr_stats_struct.facename_stats.F.pearson_std_dev',corr_stats_struct.facename_stats.M.pearson_std_dev');

    fig2 = figure();
    set(gcf, 'Position', fig_position_size,'color','w'); % fullforppt size
    t = tiledlayout(1,3,'TileSpacing','compact');
    
    % facename subplot
    nexttile
    data = facename_FM_arr;
    err = facename_FM_err_arr;

    b4 = bar(data, 'facecolor', 'flat');
    b4(1).FaceColor = 'r';
    b4(2).FaceColor = 'b';
    ax4 = gca;
    set(ax4,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    if ~permtest_or_nah
        ylim(ylims);
    end
    title('Pearson R of Median-Performing FN-TR Model','FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    [ngroups,nbars] = size(data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b4(i).XEndPoints;
    end
    
    hold on
    er = errorbar(x',data,err,'k','linestyle','none');
    hold off
    
    % ravlt_L subplot
    nexttile
    data = ravlt_L_FM_arr;
    err = ravlt_L_FM_err_arr;

    b1 = bar(data, 'facecolor', 'flat'); 
    b1(1).FaceColor = 'r';
    b1(2).FaceColor = 'b';
    ax1 = gca; % gca = get current axes
    set(ax1,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title('Pearson R of Median-Performing RAVLT-L Model','FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    [ngroups,nbars] = size(data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b1(i).XEndPoints;
    end
    
    hold on
    er = errorbar(x',data,err,'k','linestyle','none');
    hold off
    
    % ravlt_IR subplot
    nexttile
    data = ravlt_IR_FM_arr;
    err = ravlt_IR_FM_err_arr;

    b2 = bar(data, 'facecolor', 'flat'); 
    b2(1).FaceColor = 'r';
    b2(2).FaceColor = 'b';
    ax2 = gca; % gca = get current axes
    set(ax2,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title('Pearson R of Median-Performing RAVLT-IR Model','FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    [ngroups,nbars] = size(data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b2(i).XEndPoints;
    end
    
    hold on
    er = errorbar(x',data,err,'k','linestyle','none');
    hold off
    
    % link y-axes of all subplots
    linkaxes([ax4 ax1 ax2], 'y');
    
    % set custom legend for each sex
    lgd = legend(b2, {'F','M'},'Location', 'eastoutside');

    % save figure
    if ~permtest_or_nah
        saveas(fig2,'../cpm_figures/R_comparison_figs/median/median_pearson_comparison_by_sex_wErrorbars.png')
    end
    
    if permtest_or_nah
        saveas(fig2,'../cpm_figures/R_comparison_figs/median/median_pearson_comparison_by_sex_permtest_wErrorbars.png')
    end
end

if strcmp(stats_type, 'spearmans')
    %% whole group 
    ylims = [-0.05 0.5];

    fig = figure();
    set(gcf, 'Position', fig_position_size,'color','w'); % fullforpptsize
    t = tiledlayout(1,3,'TileSpacing','compact');
    
    % facename subplot
    nexttile
    data = corr_stats_struct.facename_stats.all.spearmans_median_coeff';
    errlow = corr_stats_struct.facename_stats.all.spearmans_std_dev';
    errhigh = corr_stats_struct.facename_stats.all.spearmans_std_dev';

    b4 = bar(data, 'facecolor', [0.4940 0.1840 0.5560]);
    ax4 = gca;
    set(ax4,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    if ~permtest_or_nah
        ylim(ylims);
    end
    title("Spearman's rho of Median-Performing FN-TR Model",'FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    hold on
    er = errorbar(1:numel(scan_type_list),data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    
    % ravlt_L subplot
    nexttile
    data = corr_stats_struct.ravlt_L_stats.all.spearmans_median_coeff';
    errlow = corr_stats_struct.ravlt_L_stats.all.spearmans_std_dev';
    errhigh = corr_stats_struct.ravlt_L_stats.all.spearmans_std_dev';

    b1 = bar(data, 'facecolor', [0.4940 0.1840 0.5560]);
    ax1 = gca;
    set(ax1,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title("Spearman's rho of Median-Performing RAVLT-L Model",'FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    hold on
    er = errorbar(1:numel(scan_type_list),data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    
    % ravlt_IR subplot
    nexttile
    data = corr_stats_struct.ravlt_IR_stats.all.spearmans_median_coeff';
    errlow = corr_stats_struct.ravlt_IR_stats.all.spearmans_std_dev';
    errhigh = corr_stats_struct.ravlt_IR_stats.all.spearmans_std_dev';

    b2 = bar(data, 'facecolor', [0.4940 0.1840 0.5560]);
    ax2 = gca;
    set(ax2,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title("Spearman's rho of Median-Performing RAVLT-IR Model",'FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    hold on
    er = errorbar(1:numel(scan_type_list),data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    
    % link y-axes of all subplots
    linkaxes([ax4 ax1 ax2], 'y');
    
    % save figure
    if ~permtest_or_nah
        saveas(fig,'../cpm_figures/R_comparison_figs/median/median_spearmans_comparison_all_subjs_wErrorbars.png')
    end
    
    if permtest_or_nah
        saveas(fig,'../cpm_figures/R_comparison_figs/median/median_spearmans_comparison_all_subjs_permtest_wErrorbars.png')
    end
    
    %% sex-based group 
    ylims = [-0.06 0.55];
    
    % combine female and male arrays to create stacked bar graphs
    ravlt_L_FM_arr = cat(2, corr_stats_struct.ravlt_L_stats.F.spearmans_median_coeff',corr_stats_struct.ravlt_L_stats.M.spearmans_median_coeff');
    ravlt_IR_FM_arr = cat(2, corr_stats_struct.ravlt_IR_stats.F.spearmans_median_coeff',corr_stats_struct.ravlt_IR_stats.M.spearmans_median_coeff');
    facename_FM_arr = cat(2, corr_stats_struct.facename_stats.F.spearmans_median_coeff',corr_stats_struct.facename_stats.M.spearmans_median_coeff');

    ravlt_L_FM_err_arr = cat(2, corr_stats_struct.ravlt_L_stats.F.spearmans_std_dev',corr_stats_struct.ravlt_L_stats.M.spearmans_std_dev');
    ravlt_IR_FM_err_arr = cat(2, corr_stats_struct.ravlt_IR_stats.F.spearmans_std_dev',corr_stats_struct.ravlt_IR_stats.M.spearmans_std_dev');
    facename_FM_err_arr = cat(2, corr_stats_struct.facename_stats.F.spearmans_std_dev',corr_stats_struct.facename_stats.M.spearmans_std_dev');
    
    fig2 = figure();
    set(gcf, 'Position', fig_position_size,'color','w'); % fullforppt size
    t = tiledlayout(1,3,'TileSpacing','compact');
    
    % facename subplot
    nexttile
    data = facename_FM_arr;
    err = facename_FM_err_arr;

    b4 = bar(data, 'facecolor', 'flat');
    b4(1).FaceColor = 'r';
    b4(2).FaceColor = 'b';
    ax4 = gca;
    set(ax4,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    if ~permtest_or_nah
        ylim(ylims);
    end
    title("Spearman's rho of Median-Performing FN-TR Model",'FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);
    
    [ngroups,nbars] = size(data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b4(i).XEndPoints;
    end
    
    hold on
    er = errorbar(x',data,err,'k','linestyle','none');
    hold off
    
    % ravlt_L subplot
    nexttile
    data = ravlt_L_FM_arr;
    err = ravlt_L_FM_err_arr;

    b1 = bar(data, 'facecolor', 'flat'); 
    b1(1).FaceColor = 'r';
    b1(2).FaceColor = 'b';
    ax1 = gca; % gca = get current axes
    set(ax1,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title("Spearman's rho of Median-Performing RAVLT-L Model",'FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);
    
    [ngroups,nbars] = size(data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b1(i).XEndPoints;
    end
    
    hold on
    er = errorbar(x',data,err,'k','linestyle','none');
    hold off
    
    % ravlt_IR subplot
    nexttile
    data = ravlt_IR_FM_arr;
    err = ravlt_IR_FM_err_arr;

    b2 = bar(data, 'facecolor', 'flat'); 
    b2(1).FaceColor = 'r';
    b2(2).FaceColor = 'b';
    ax2 = gca; % gca = get current axes
    set(ax2,'XTick',1:numel(scan_type_list),'XTickLabel',scan_type_list,'TickLabelInterpreter', 'none','FontName','Times New Roman','fontsize',fs);
    title("Spearman's rho of Median-Performing RAVLT-IR Model",'FontName','Times New Roman','fontsize',title_fs);
    xtickangle(30);

    [ngroups,nbars] = size(data);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b2(i).XEndPoints;
    end
    
    hold on
    er = errorbar(x',data,err,'k','linestyle','none');
    hold off
    
    % link y-axes of all subplots
    linkaxes([ax4 ax1 ax2], 'y');
    
    % set custom legend for each sex
    lgd = legend(b2, {'F','M'},'Location', 'eastoutside');

    % save figure
    if ~permtest_or_nah
        saveas(fig2,'../cpm_figures/R_comparison_figs/median/median_spearmans_comparison_by_sex_wErrorbars.png')
    end
    
    if permtest_or_nah
        saveas(fig2,'../cpm_figures/R_comparison_figs/median/median_spearmans_comparison_by_sex_permtest_wErrorbars.png')
    end
end

end

%% fxn to perform Wilcoxon rank sum test
function [p_val, z_val, rs] = wilcoxon_rs_test(F_stats_arr, M_stats_arr, n_runs)

[p,h,stats] = ranksum(F_stats_arr(1,:),M_stats_arr(1,:));
p_val = p;
z_val = stats.zval;
rs = stats.ranksum;

end