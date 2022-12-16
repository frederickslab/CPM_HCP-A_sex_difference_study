%% Written by Suyeon Ju, 7.10.22, adapted from Corey Horien's scripts
%% Implementation
function pmask_visualization(pos_mat,neg_mat,param_name,scan_type_name, group, cov,standardized_or_nah)
%% visualization setup
no_nodes = 268;
no_networks = 10;

ten_network_defn_path =  '../shen_268_labels';
filename = 'ten_network_defn.mat';
file = fullfile(ten_network_defn_path, filename);
load(file);

network_ticklabels = {'MF','FP','DMN','Mot','VI','VII','VAs','SAL','SC','CBL'};

% *** CHECK FILENAMES (ESP K AND ITER THRESHOLDS) BEFORE RUNNING!!!

% adjust colorbar value limits to standardized plots within each param
if strcmp(param_name, 'facename')
    pos_mat_y_lims_standardized = [0 0.185];
    neg_mat_y_lims_standardized = [0 0.15];
elseif strcmp(param_name, 'ravlt_L')
    pos_mat_y_lims_standardized = [0 0.14];
    neg_mat_y_lims_standardized = [0 0.125];
elseif strcmp(param_name, 'ravlt_IR')
    pos_mat_y_lims_standardized = [0 0.17];
    neg_mat_y_lims_standardized = [0 0.185];
end

%% loading in files you want to group by network size. <- just use pos_mat + neg_mat for now!

% Now do whatever you want with this file name, such as reading it in as an image array with imread()
mats = cat(3,pos_mat,neg_mat);

for k = 1 : size(mats,3)
    thr_1 = mats(:,:,k);
    %this is arranging the networks in the original matrix (into ten diff networks), first by rows, then by columns.
    new_assignments_matr_INTERMEDIATE = thr_1(ten_network_defn(:,2),:);
    new_assignments_final_matr_tmp = new_assignments_matr_INTERMEDIATE(:,ten_network_defn(:,2));
    new_assignments_final_matr_cell{k} = new_assignments_final_matr_tmp;
end

%% this is finding the number of sig edges across the entire 10x10 network


for w = 1:length(new_assignments_final_matr_cell);
    
    new_assignments_final_matr = new_assignments_final_matr_cell{1,w};
    
       
    for mm = 1:no_networks
        for k = 1:no_networks
            zero_matrix = zeros(no_nodes, no_nodes);
            [indices] = find( ten_network_defn(:, 3)==k); % this will indicate the row to be used in the between network comparison
            indices_network_mm = find(ten_network_defn(:,3) ==mm); %this is indicating the column to be used in the between network comparison

            zero_matrix(indices,indices_network_mm) = 1;
            network_k = zero_matrix;
            number_of_edges = (new_assignments_final_matr + network_k);

            tmp_raw_DP_edges_within_network_mm = length(find(number_of_edges == 2));


            %GETTING RAW NUMBER OF EDGES
            if mm ==k
                raw_DP_edges_within_network_mm(k) = tmp_raw_DP_edges_within_network_mm/(2);
            elseif mm ~=k
                raw_DP_edges_within_network_mm(k) = tmp_raw_DP_edges_within_network_mm;
            end


            [indices_size,~] = size(indices); %getting the size of the network so I can use it to normalize below. --> this is what I'll need to change for SLIM and UPSM.
            [indices_network_mm_size,~] = size(indices_network_mm); %getting the size of the network so I can use it to normalize below.

            %GETTING EDGES/NETWORK SIZE

            if mm == k %note that this is for the case when I am dividing by the within network edges.
                edges_divided_by_net_size(k) =  (tmp_raw_DP_edges_within_network_mm/(2))/((indices_size*indices_network_mm_size - indices_size)./2); %here I'm just dividing the number of edges in that square by the total network size.
            elseif mm ~= k %note that this is for the case when I am dividing by the between-network edges.
                edges_divided_by_net_size(k) =  tmp_raw_DP_edges_within_network_mm./(indices_size*indices_network_mm_size);
            end

        end

        mat_test_raw_edges(mm,:) = raw_DP_edges_within_network_mm; %RAW NUMBER OF EDGES
        mat_test_edges_by_net_size(mm,:) = edges_divided_by_net_size; %GETTING EDGES/NETWORK SIZE

    end


    mat_test_raw_edges_lower_tri = tril(mat_test_raw_edges,0);
    mat_test_edges_by_net_size_lower_tri =  tril(mat_test_edges_by_net_size,0);

    mat_1{w} = mat_test_raw_edges_lower_tri; % mat_1 is lower triangle of raw edges
    mat_2{w} = mat_test_edges_by_net_size_lower_tri; % mat_2 is lower triangle of edges normalized by network size
    
end


%% note that these figures are for quick visualization only. they should be
%tinkered  with or a different script should be used for
%publication/presentation purposes.

% %mat_1 --> plotting the number of edges in a network without normalizing by
% %network size.
% figure;
% subplot(2,3,1)
% imagesc(mat_1{1,1})
% colorbar
% title('pos mask')
% 
% subplot(2,3,2)
% imagesc(mat_1{1,2})
% colorbar
% title('neg mask')
% 
% % raw_edge_fig_filename = sprintf('/Users/sj737/Library/CloudStorage/OneDrive-YaleUniversity/Fredericks_Lab_files/CPM_HCP-A/CPM_HCP-Aging/figs_and_csvmats/p%.2f_k%d/%s/consensus_mask_figs/pos_neg_mask_network_representation_raw_edges_%s_%s_p%.2f_k%d.png', p_thresh, k_folds, param, scan_type, param, p_thresh, k_folds);
% % raw_edge_fig_filename = sprintf('/Users/sj737/Library/CloudStorage/OneDrive-YaleUniversity/Fredericks_Lab_files/CPM_HCP-A/CPM_HCP-Aging/figs_and_csvmats/p%.3f_k%d/%s/consensus_mask_figs/pos_neg_mask_network_representation_raw_edges_%s_%s_p%.3f_k%d.png', p_thresh, k_folds, param, scan_type, param, p_thresh, k_folds);
% % saveas(gcf,raw_edge_fig_filename)
% 
% 
% %mat_2 --> plotting the number of edges in a network normalizing by
% %network size.
% 
% figure;
% subplot(2,3,1)
% imagesc(mat_2{1,1})
% colorbar
% title('pos mask')
% 
% subplot(2,3,2)
% imagesc(mat_2{1,2})
% colorbar
% title('neg mask')
% 
% % norm_edge_fig_filename = sprintf('/Users/sj737/Library/CloudStorage/OneDrive-YaleUniversity/Fredericks_Lab_files/CPM_HCP-A/CPM_HCP-Aging/figs_and_csvmats/p%.2f_k%d/%s/consensus_mask_figs/pos_neg_mask_network_representation_raw_edges_normalized_%s_%s_p%.2f_k%d.png', p_thresh, k_folds, param, scan_type, param, p_thresh, k_folds);
% % norm_edge_fig_filename = sprintf('/Users/sj737/Library/CloudStorage/OneDrive-YaleUniversity/Fredericks_Lab_files/CPM_HCP-A/CPM_HCP-Aging/figs_and_csvmats/p%.3f_k%d/%s/consensus_mask_figs/pos_neg_mask_network_representation_raw_edges_normalized_%s_%s_p%.3f_k%d.png', p_thresh, k_folds, param, scan_type, param, p_thresh, k_folds);
% % saveas(gcf,norm_edge_fig_filename)

%% mat_2_1 - plotting pos mat of mat_2 (lower triangle of edges normalized by network size)

SLIM_DP_edges = mat_2{1,1};

SLIM_DP_edges(:,end+1) = 0; % pad with zeros for pcolor
SLIM_DP_edges(end+1,:) = 0;

%finding the max value of each matrix so I can use it to scale the color bar below.
tmp_max = max(SLIM_DP_edges);
max_SLIM = max(tmp_max);
% ^I'm putting the max value I obtained from directly above in the zero spaces
% of the upper triangle in my matrix -- I am doing this because they are
% currently filled with zeros, which results in that square becoming black
% in my final figure with the current colormap. 

filler_upper_triangle_SLIM = zeros(10,10)*max_SLIM;
filler_upper_triangle_SLIM = triu(filler_upper_triangle_SLIM,1);

filler_upper_triangle_SLIM =[filler_upper_triangle_SLIM zeros(10,1)]; %here I am adding zeros so I can add filler_upper_triangle_SLIM to my matrices below.
filler_upper_triangle_SLIM = vertcat(filler_upper_triangle_SLIM, zeros(1,11)); 

%putting them together
final_SLIM = SLIM_DP_edges + filler_upper_triangle_SLIM;

% Plot and format figure
figure;
pcolor(final_SLIM);
if ~standardized_or_nah
    caxis([0 max_SLIM]); %color bar will be scaled from 0 to the max value in the matrix.
end
if standardized_or_nah
    caxis(pos_mat_y_lims_standardized);
end
axis('square');
set(gca,'YDir','reverse','XTickLabel',network_ticklabels,'YTickLabel',network_ticklabels,'Xtick',1.5 : 1 : 10.5,'Ytick',1.5 : 1 : 10.5);
set(gcf,'color','w');
colormap(brewermap(256,'reds'))
% colormap hot
colorbar
% title(sprintf("%s - %s - %s positive edges", param_name, scan_type_name, group), 'Interpreter', 'none');
% title(sprintf("%s - %s 10-network positive-edge heatmap", param_name, scan_type_name), 'Interpreter', 'none');
set(gcf, 'Position', [1000, 500, 300, 250],'color','w');

if ~standardized_or_nah
    pos_mat_vis_fig_filename = sprintf('../cpm_figures/ten_network_heatmaps/use_for_paper/%s_%s_pos_mat_heatmap_k40_iter40_%s_%s.png', param_name, scan_type_name, group, cov);
end
if standardized_or_nah
    pos_mat_vis_fig_filename = sprintf('../cpm_figures/ten_network_heatmaps/use_for_paper/%s_%s_pos_mat_heatmap_k40_iter40_%s_%s_standardized.png', param_name, scan_type_name, group, cov);
end

saveas(gcf,pos_mat_vis_fig_filename)

%% mat_2_2 - plotting neg mat of mat_2 (lower triangle of edges normalized by network size)
SLIM_DP_edges = mat_2{1,2};

SLIM_DP_edges(:,end+1) = 0; % pad with zeros for pcolor
SLIM_DP_edges(end+1,:) = 0;

%finding the max value of each matrix so I can use it to scale the color
%bar below.

tmp_max = max(SLIM_DP_edges);
max_SLIM = max(tmp_max);
% max_SLIM = 0;

%I'm putting the max value I obtained from directly above in the zero spaces
%of the upper triangle in my matrix -- I am doing this because they are
%currently filled with zeros, which results in that square becoming black
%in my final figure with the current colormap. 

filler_upper_triangle_SLIM = zeros(10,10)*max_SLIM;
filler_upper_triangle_SLIM = triu(filler_upper_triangle_SLIM,1);

filler_upper_triangle_SLIM =[filler_upper_triangle_SLIM zeros(10,1)]; %here I am adding zeros so I can add filler_upper_triangle_SLIM to my matrices below.
filler_upper_triangle_SLIM = vertcat(filler_upper_triangle_SLIM, zeros(1,11)); 

%putting them together
final_SLIM = SLIM_DP_edges + filler_upper_triangle_SLIM;

% Plot and format figure
figure;
pcolor(final_SLIM);
if ~standardized_or_nah
    caxis([0 max_SLIM]); %color bar will be scaled from 0 to the max value in the matrix.
end
if standardized_or_nah
    caxis(neg_mat_y_lims_standardized);
end
axis('square');
set(gca,'YDir','reverse','XTickLabel',network_ticklabels,'YTickLabel',network_ticklabels,'Xtick',1.5 : 1 : 10.5,'Ytick',1.5 : 1 : 10.5);
set(gcf,'color','w');
colormap(brewermap(256,'blues'))
% colormap hot
colorbar
% title(sprintf("%s - %s - %s negative edges", param_name, scan_type_name, group), 'Interpreter', 'none');
% title(sprintf("%s - %s 10-network negative-edge heatmap", param_name, scan_type_name), 'Interpreter', 'none');
set(gcf, 'Position', [1000, 100, 300, 250],'color','w');

if ~standardized_or_nah
    neg_mat_vis_fig_filename = sprintf('../cpm_figures/ten_network_heatmaps/use_for_paper/%s_%s_neg_mat_heatmap_k40_iter40_%s_%s.png', param_name, scan_type_name, group, cov);
end
if standardized_or_nah
    neg_mat_vis_fig_filename = sprintf('../cpm_figures/ten_network_heatmaps/use_for_paper/%s_%s_neg_mat_heatmap_k40_iter40_%s_%s_standardized.png', param_name, scan_type_name, group, cov);
end
saveas(gcf,neg_mat_vis_fig_filename)

clear all

end