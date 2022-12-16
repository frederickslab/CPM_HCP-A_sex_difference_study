%% Written by Suyeon Ju, 7.10.22, adapted from Corey Horien's scripts

%% general script info

% fxn to get positive and negative pmasks of significant edges from cpm

% inputs:
%   `param_list` = cell array of parameters to be tested
%       i.e., "{'ravlt','neon'}"

% outputs:
%   `MxM_matrix_pos`
%   `MxM_matrix_neg`
%   `size_of_pos_mask`
%   `size_of_neg_mask`

%% Implementation
function [MxM_matrix_pos,MxM_matrix_neg,size_of_pos_mask,size_of_neg_mask] = get_consensus_mask(pmask_stock,k_folds,trial_count,k_fold_threshold,iter_threshold)

k_thresh = k_fold_threshold/k_folds;
iter_thresh = iter_threshold/trial_count;

%% condense pmask_stock (across all folds and all trials)
% average pmask values across all folds (add up pmask values, then divide by k_folds number)
pmask_fold = (sum(pmask_stock, 2))/k_folds;

% set all edges with magnitude less than k_thresh (aka, less than 3 out of 5 folds = 0.6) to zero, and all others to 1 or -1 (depending on sign of average)
pmask_fold(pmask_fold < k_thresh & pmask_fold > -k_thresh) = 0;
pmask_fold(pmask_fold >= k_thresh) = 1;
pmask_fold(pmask_fold <= (-1*k_thresh)) = -1;

% average pmask values across all trials
pmask_trial = sum(squeeze(pmask_fold),2)/trial_count;

% set all edges with magnitude less than iter_thresh (aka, less than 600 out of 1000 iterations = 0.6) to zero, and all others to 1 or -1 (depending on sign of average)
pmask_trial(pmask_trial < iter_thresh & pmask_trial > -iter_thresh) = 0;
pmask_trial(pmask_trial >= iter_thresh) = 1;
pmask_trial(pmask_trial <= (-1*iter_thresh)) = -1;

pmask = pmask_trial;

% extract positive and negative elements from pmask
pmask_pos = pmask;
pmask_pos(pmask_pos == -1) = 0;
% pmask_pos(pmask_pos<=0) = 0;

pmask_neg = pmask;
pmask_neg(pmask_neg == 1) = 0;
% pmask_neg(pmask_neg>=0) = 0;

%% check checking that there are no values in the mats that shouldn't be there
find(pmask_pos == -1);
find(pmask_neg == 1);

%% these are the sizes of edges that show up 
size_of_pos_mask = length(find(pmask_pos));
size_of_neg_mask = length(find(pmask_neg));

%% triangularization of pos_mat/neg_mat vectors
no_node = 268;
aa = ones(no_node, no_node);
aa_upp = triu(aa, 1);
upp_id = find(aa_upp);
upp_len = length(upp_id);

%back in matrix format
edge_vector_matrix_pos = zeros(no_node, no_node);
edge_vector_matrix_neg = zeros(no_node, no_node);

%now need to put idx_pos into 35,778 vector - makes the 268x268 pmask
%without removing any node data (will do this later with thresholding)
edge_vector_matrix_pos(upp_id) = pmask_pos;
edge_vector_matrix_pos = edge_vector_matrix_pos + edge_vector_matrix_pos';
MxM_matrix_pos = edge_vector_matrix_pos;

edge_vector_matrix_neg(upp_id) = pmask_neg;
edge_vector_matrix_neg = edge_vector_matrix_neg + edge_vector_matrix_neg';
MxM_matrix_neg = edge_vector_matrix_neg * -1;
clear pmask_pos pmask_neg

end