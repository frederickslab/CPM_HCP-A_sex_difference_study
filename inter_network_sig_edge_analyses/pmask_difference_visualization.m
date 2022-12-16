%% Written by Suyeon Ju, 7.10.22, adapted from Corey Horien's scripts
%% Implementation
function pmask_difference_visualization(F_mat,M_mat,param_name,scan_type_name,cov,standardized_or_nah, pos_or_neg)
%% visualization setup
no_nodes = 268;
no_networks = 10;

ten_network_defn_path = '../shen_268_labels';
filename = 'ten_network_defn.mat';
file = fullfile(ten_network_defn_path, filename);
load(file);

network_ticklabels = {'MF','FP','DMN','Mot','VI','VII','VAs','SAL','SC','CBL'};

FM_mats = cat(3,F_mat,M_mat);

for k = 1 : size(FM_mats,3)
    thr_1 = FM_mats(:,:,k);
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
            [indices] = find( ten_network_defn(:, 3) == k); % this will indicate the row to be used in the between network comparison
            indices_network_mm = find(ten_network_defn(:,3) == mm); %this is indicating the column to be used in the between network comparison

            zero_matrix(indices,indices_network_mm) = 1;
            network_k = zero_matrix;
            number_of_edges = (new_assignments_final_matr + network_k);

            tmp_raw_DP_edges_within_network_mm = length(find(number_of_edges == 2));

            %GETTING RAW NUMBER OF EDGES
            if mm == k
                raw_DP_edges_within_network_mm(k) = tmp_raw_DP_edges_within_network_mm/(2);
            elseif mm ~= k
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

%% plot F-M difference pmask (lower triangle of edges normalized by network size)

FM_difference_ten_network_edges = mat_2{1,1}-mat_2{1,2};
% size(FM_difference_ten_network_edges)

FM_difference_ten_network_edges(:,end+1) = 0; % pad with zeros for pcolor
FM_difference_ten_network_edges(end+1,:) = 0;

%finding the max value of each matrix so I can use it to scale the color bar below.
tmp_max = max(FM_difference_ten_network_edges);
max_FM_diff = max(tmp_max);

tmp_min = min(FM_difference_ten_network_edges);
min_FM_diff = min(tmp_min);

if max_FM_diff > abs(min_FM_diff)
    color_lim = max_FM_diff;
else
    color_lim = abs(min_FM_diff);
end

filler_upper_triangle_FM_diff = zeros(10,10)*max_FM_diff;
filler_upper_triangle_FM_diff = triu(filler_upper_triangle_FM_diff,1);

filler_upper_triangle_FM_diff =[filler_upper_triangle_FM_diff zeros(10,1)]; %here I am adding zeros so I can add filler_upper_triangle_SLIM to my matrices below.
filler_upper_triangle_FM_diff = vertcat(filler_upper_triangle_FM_diff, zeros(1,11)); 

%putting them together
final_FM_diff = FM_difference_ten_network_edges + filler_upper_triangle_FM_diff;

final_FM_diff(final_FM_diff == 0) = NaN;

% Plot and format figure
figure;
pcolor(final_FM_diff);
if ~standardized_or_nah
    caxis([-color_lim color_lim]); %color bar will be scaled from 0 to the max value in the matrix.
else
    caxis(mat_y_lims_standardized);
end
axis('square');
set(gca,'YDir','reverse','XTickLabel',network_ticklabels,'YTickLabel',network_ticklabels,'Xtick',1.5 : 1 : 10.5,'Ytick',1.5 : 1 : 10.5);
set(gcf,'color','w');
colormap(brewermap(256,'-RdBu'))
% colormap hot
colorbar
% title(sprintf("%s - %s F-M difference 10-network %s edges", param_name, scan_type_name, pos_or_neg), 'Interpreter', 'none');
set(gcf, 'Position', [1000, 500, 300, 250],'color','w');

if ~standardized_or_nah
    mat_vis_fig_filename = sprintf('../cpm_figures/ten_network_heatmaps/use_for_paper/%s_%s_heatmap_k40_iter40_%s_%s_mat.png', param_name, scan_type_name, cov, pos_or_neg);
end
if standardized_or_nah
    mat_vis_fig_filename = sprintf('../cpm_figures/ten_network_heatmaps/use_for_paper/%s_%s_heatmap_k40_iter40_%s_%s_mat_standardized.png', param_name, scan_type_name, cov, pos_or_neg);
end

saveas(gcf,mat_vis_fig_filename)

end