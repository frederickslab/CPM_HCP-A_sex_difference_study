% intra-DMN pmask visualization script

close all
clear all

%% setup variables
% adjust the variables below as needed!
param_list = {'facename','ravlt_L','ravlt_IR'};
scan_type = 'tfMRI_FACENAME';
standardized_or_nah = 1; % 0 = colorbars NOT standardized between all figs; 1 = colorbars ARE standardized between all figs
plot_F_and_M_pmasks_separately_bool = 0; % aka, NOT plotting the difference of F and M pmasks!
plot_F_and_M_pmasks_difference_bool = 1;
plot_intra_DMN_summed_vec_bool = 0;

%% plot all pmasks as heatmaps
for i = 1:length(param_list)
    param = param_list{i};
    % load in DMN pmasks for each behav parameter (M x M x 4 matrix for each param, where M is 20 nodes, and the 4 are listed below)
    %   DMN_matrix(:,:,1) = pos intra-DMN pmask of M model
    %   DMN_matrix(:,:,2) = neg intra-DMN pmask of M model
    %   DMN_matrix(:,:,3) = pos intra-DMN pmask of F model
    %   DMN_matrix(:,:,4) = neg intra-DMN pmask of F model
    
    if plot_F_and_M_pmasks_separately_bool
        load(sprintf('../BIG_data_from_CPM_HCP-Aging/intra-DMN-pmasks/%s_DMNedges_sorted.mat', param_list{i}))
        % plot pos intra-DMN pmask of M model
        plot_intra_DMN_pos_mat(DMN_matrix(:,:,1), param, scan_type, 'M', standardized_or_nah);
        % plot neg intra-DMN pmask of M model
        plot_intra_DMN_neg_mat(-DMN_matrix(:,:,2), param, scan_type, 'M', standardized_or_nah);
    
        % plot pos intra-DMN pmask of F model
        plot_intra_DMN_pos_mat(DMN_matrix(:,:,3), param, scan_type, 'F', standardized_or_nah);
        % plot neg intra-DMN pmask of F model
        plot_intra_DMN_neg_mat(-DMN_matrix(:,:,4), param, scan_type, 'F', standardized_or_nah);
    end

    if plot_F_and_M_pmasks_difference_bool
        load(sprintf('../BIG_data_from_CPM_HCP-Aging/intra-DMN-pmasks/%s_DMNedges_sorted.mat', param_list{i}))
        % plot F-M difference pos intra-DMN pmask
        pos_pmask_F_M_difference = DMN_matrix(:,:,3) - DMN_matrix(:,:,1);
        plot_intra_DMN_difference(pos_pmask_F_M_difference, param, scan_type, 'positive', standardized_or_nah);

        % plot F-M difference neg intra-DMN pmask
        neg_pmask_F_M_difference = DMN_matrix(:,:,4) - DMN_matrix(:,:,2);
        plot_intra_DMN_difference(neg_pmask_F_M_difference, param, scan_type, 'negative', standardized_or_nah);
    end

    if plot_intra_DMN_summed_vec_bool
        load(sprintf('../BIG_data_from_CPM_HCP-Aging/intra-DMN-pmasks/%s_DMNedges_summed.mat', param_list{i}))
        % plot F-M difference pos intra-DMN summed node vector
        pos_sum_vec_F_M_difference = DMN_summed(:,3) - DMN_summed(:,1);
        plot_intra_DMN_summed_vec(pos_sum_vec_F_M_difference, param, scan_type, 'positive', standardized_or_nah);

        % plot F-M difference neg intra-DMN summed node vector
        neg_sum_vec_F_M_difference = DMN_summed(:,4) - DMN_summed(:,2);
        plot_intra_DMN_summed_vec(neg_sum_vec_F_M_difference, param, scan_type, 'negative', standardized_or_nah);
    end

end


%% fxn to plot pos mat of model (lower triangle of edges normalized by network size)
function plot_intra_DMN_pos_mat(dmn_mat, param, scan_type, group, standardized_or_nah)

% adjust colorbar value limits to standardized plots within each param
if strcmp(param, 'facename')
    pos_mat_y_lims_standardized = [0 5000]; 
elseif strcmp(param, 'ravlt_L')
    pos_mat_y_lims_standardized = [0 4900];
elseif strcmp(param, 'ravlt_IR')
    pos_mat_y_lims_standardized = [0 5000];
end

intra_DMN_edges = dmn_mat;

intra_DMN_edges(:,end+1) = 0; % pad with zeros for pcolor
intra_DMN_edges(end+1,:) = 0;

%finding the max value of each matrix so I can use it to scale the color bar below.
tmp_max = max(intra_DMN_edges);
max_intra_DMN = max(tmp_max);
% ^I'm putting the max value I obtained from directly above in the zero spaces
% of the upper triangle in my matrix -- I am doing this because they are
% currently filled with zeros, which results in that square becoming black
% in my final figure with the current colormap. 

filler_upper_triangle_intra_DMN = zeros(20,20)*max_intra_DMN;
filler_upper_triangle_intra_DMN = triu(filler_upper_triangle_intra_DMN,1);

filler_upper_triangle_intra_DMN =[filler_upper_triangle_intra_DMN zeros(20,1)]; %here I am adding zeros so I can add filler_upper_triangle_SLIM to my matrices below.
filler_upper_triangle_intra_DMN = vertcat(filler_upper_triangle_intra_DMN, zeros(1,21)); 

%putting them together
final_intra_DMN = intra_DMN_edges + filler_upper_triangle_intra_DMN;
final_intra_DMN = tril(final_intra_DMN,0);

% Plot and format figure
figure;
pcolor(final_intra_DMN);
if ~standardized_or_nah
    caxis([0 max_intra_DMN]); %color bar will be scaled from 0 to the max value in the matrix.
else
    caxis(pos_mat_y_lims_standardized);
end
axis('square');
set(gca,'YDir','reverse','XTickLabel',1:1:20,'YTickLabel',1:1:20,'Xtick',1.5 : 1 : 20.5,'Ytick',1.5 : 1 : 20.5);
set(gcf,'color','w');
colormap(brewermap(256,'reds'))
% colormap hot
colorbar
title(sprintf("%s - %s %s-group intra-DMN positive edges", param, scan_type, group), 'Interpreter', 'none');
set(gcf, 'Position', [1000, 500, 500, 500],'color','w');

if ~standardized_or_nah
    pos_mat_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/%s_pos_mat_%s_%s_intra_DMN_heatmap.png', param, scan_type, group);
else
    pos_mat_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/%s_pos_mat_%s_%s_intra_DMN_heatmap_standardized.png', param, scan_type, group);
end

saveas(gcf,pos_mat_vis_fig_filename)
end

%% fxn to plot pos mat of model (lower triangle of edges normalized by network size)
function plot_intra_DMN_neg_mat(dmn_mat, param, scan_type, group, standardized_or_nah)

% adjust colorbar value limits to standardized plots within each param
if strcmp(param, 'facename')
    neg_mat_y_lims_standardized = [0 2750];
elseif strcmp(param, 'ravlt_L')
    neg_mat_y_lims_standardized = [0 2400];
elseif strcmp(param, 'ravlt_IR')
    neg_mat_y_lims_standardized = [0 2350];
end

intra_DMN_edges = dmn_mat;

intra_DMN_edges(:,end+1) = 0; % pad with zeros for pcolor
intra_DMN_edges(end+1,:) = 0;

%finding the max value of each matrix so I can use it to scale the color bar below.
tmp_max = max(intra_DMN_edges);
max_intra_DMN = max(tmp_max);
% ^I'm putting the max value I obtained from directly above in the zero spaces
% of the upper triangle in my matrix -- I am doing this because they are
% currently filled with zeros, which results in that square becoming black
% in my final figure with the current colormap. 

filler_upper_triangle_intra_DMN = zeros(20,20)*max_intra_DMN;
filler_upper_triangle_intra_DMN = triu(filler_upper_triangle_intra_DMN,1);

filler_upper_triangle_intra_DMN =[filler_upper_triangle_intra_DMN zeros(20,1)]; %here I am adding zeros so I can add filler_upper_triangle_SLIM to my matrices below.
filler_upper_triangle_intra_DMN = vertcat(filler_upper_triangle_intra_DMN, zeros(1,21)); 

%putting them together
final_intra_DMN = intra_DMN_edges + filler_upper_triangle_intra_DMN;
final_intra_DMN = tril(final_intra_DMN,0);

% Plot and format figure
figure;
pcolor(final_intra_DMN);
if ~standardized_or_nah
    caxis([0 max_intra_DMN]); %color bar will be scaled from 0 to the max value in the matrix.
else
    caxis(neg_mat_y_lims_standardized);
end
axis('square');
set(gca,'YDir','reverse','XTickLabel',1:1:20,'YTickLabel',1:1:20,'Xtick',1.5 : 1 : 20.5,'Ytick',1.5 : 1 : 20.5);
set(gcf,'color','w');
colormap(brewermap(256,'blues'))
% colormap hot
colorbar
title(sprintf("%s - %s %s-group intra-DMN negative edges", param, scan_type, group), 'Interpreter', 'none');
set(gcf, 'Position', [1000, 500, 500, 500],'color','w');

if ~standardized_or_nah
    neg_mat_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/%s_neg_mat_%s_%s_intra_DMN_heatmap.png', param, scan_type, group);
else
    neg_mat_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/%s_neg_mat_%s_%s_intra_DMN_heatmap_standardized.png', param, scan_type, group);
end

saveas(gcf,neg_mat_vis_fig_filename)
end

%% fxn to plot pos mat of model (lower triangle of edges normalized by network size)
function plot_intra_DMN_difference(dmn_mat, param, scan_type, pos_or_neg, standardized_or_nah)

% adjust colorbar value limits to standardized plots within each param
mat_y_lims_standardized = [-5000 5000];

intra_DMN_edges = dmn_mat;

intra_DMN_edges(:,end+1) = 0; % pad with zeros for pcolor
intra_DMN_edges(end+1,:) = 0;

% %finding the max value of each matrix so I can use it to scale the color bar below.
tmp_max = max(intra_DMN_edges);
max_intra_DMN = max(tmp_max);

tmp_min = min(intra_DMN_edges);
min_intra_DMN = min(tmp_min);

if max_intra_DMN > abs(min_intra_DMN)
    color_lim = max_intra_DMN;
else
    color_lim = abs(min_intra_DMN);
end


final_intra_DMN = intra_DMN_edges;

final_intra_DMN = tril(final_intra_DMN,0);

final_intra_DMN(final_intra_DMN == 0) = NaN;


% Plot and format figure
fig = figure;
pcolor(final_intra_DMN);

if ~standardized_or_nah
    caxis([-color_lim color_lim]); %color bar will be scaled from 0 to the max value in the matrix.
else
    caxis(mat_y_lims_standardized);
end
axis('square');
set(gca,'YDir','reverse','XTickLabel',1:1:20,'YTickLabel',1:1:20,'Xtick',1.5 : 1 : 20.5,'Ytick',1.5 : 1 : 20.5);
colormap(brewermap(256,'-RdBu'))
colorbar
% title(sprintf("%s - %s F-M difference intra-DMN %s edges", param, scan_type, pos_or_neg), 'Interpreter', 'none');
set(gcf, 'Position', [1000, 500, 500, 500],'color','w','InvertHardCopy', 'off')
set(gca,'color','w');
% yline(gca,[3 13 19],'LineWidth', 3);
% xline(gca,[3 13 19],'LineWidth', 3);

if ~standardized_or_nah
    mat_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/%s_%s_F-M_intra_DMN_%s_heatmap.png', param, scan_type, pos_or_neg);
else
    mat_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/%s_%s_F-M_intra_DMN_%s_heatmap_standardized_new.png', param, scan_type, pos_or_neg);
end

saveas(fig,mat_vis_fig_filename)
end

%% fxn to plot summed vector of intra-DMN edges
function plot_intra_DMN_summed_vec(dmn_vec, param, scan_type, pos_or_neg, standardized_or_nah)

% adjust colorbar value limits to standardized plots within each param
y_lims_standardized = [-12653 12653]; % based on largest color_lim (below)

intra_DMN_summed_edges = dmn_vec;

intra_DMN_summed_edges(:,end+1) = 0; % pad with zeros for pcolor
intra_DMN_summed_edges(end+1,:) = 0;

% %finding the max value of each matrix so I can use it to scale the color bar below.
tmp_max = max(intra_DMN_summed_edges);
max_intra_DMN = max(tmp_max);

tmp_min = min(intra_DMN_summed_edges);
min_intra_DMN = min(tmp_min);

if max_intra_DMN > abs(min_intra_DMN)
    color_lim = max_intra_DMN;
else
    color_lim = abs(min_intra_DMN);
end
% disp(color_lim) % largest color_lim is 1.265250000000000e+04

final_intra_DMN = intra_DMN_summed_edges;
final_intra_DMN(final_intra_DMN == 0) = NaN;

% Plot and format figure
fig = figure;
pcolor(final_intra_DMN);

if ~standardized_or_nah
    caxis([-color_lim color_lim]); %color bar will be scaled from 0 to the max value in the matrix.
else
    caxis(y_lims_standardized);
end
set(gca,'YDir','reverse','YTickLabel',1:1:20,'Ytick',1.5 : 1 : 20.5, 'XTick', []);
colormap(brewermap(256,'-RdBu'))
c = colorbar;
set(c, 'Location', 'eastoutside');
set(gcf, 'Position', [1000, 500, 100, 500],'color','w','InvertHardCopy', 'off')
set(gca,'color','w');

if ~standardized_or_nah
    vec_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/summed_vector_figs/%s_%s_F-M_intra_DMN_%s_summed_vector.png', param, scan_type, pos_or_neg);
else
    vec_vis_fig_filename = sprintf('../cpm_figures/intra_DMN_heatmaps/summed_vector_figs/%s_%s_F-M_intra_DMN_%s_summed_vector_standardized.png', param, scan_type, pos_or_neg);
end

saveas(fig,vec_vis_fig_filename)

end