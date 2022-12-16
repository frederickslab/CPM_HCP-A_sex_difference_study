
mpmask_facename = pmasks.M_pmasks.tfMRI_FACENAME;
fpmask_facename = pmasks.M_pmasks.tfMRI_FACENAME;

%sum across k-folds and iterations
mmatsum_kfold = sum(mpmask_facename,2);
mmatsum_tot = sum(mmatsum_kfold, 3);
fmatsum_kfold = sum(fpmask_facename ,2);
fmatsum_tot = sum(fmatsum_kfold,3);

%triangularization of mats vectors
aa = ones(268, 268);
aa_upp = triu(aa, 1);
upp_idx = find(aa_upp);
upp_len = length(upp_idx);

medge_vector_matrix = zeros(268, 268);
fedge_vector_matrix = zeros(268, 268);

medge_vector_matrix(upp_idx) = mmatsum_tot;
medge_vector_matrix = medge_vector_matrix + medge_vector_matrix';
mMxM_matrix = medge_vector_matrix;

fedge_vector_matrix(upp_idx) = mmatsum_tot;
fedge_vector_matrix = fedge_vector_matrix + fedge_vector_matrix';
fMxM_matrix = fedge_vector_matrix;

%3d matrix to put m+, m-, f+, f- 
combined_matrix = zeros(268, 268, 4);

%%
% extract positive and negative elements from pmask
mpmask_pos = mMxM_matrix;
% pmask_pos(pmask_pos == -1) = 0;
mpmask_pos(mpmask_pos<=0) = 0;
combined_matrix(:,:,1) = mpmask_pos;

mpmask_neg = mMxM_matrix;
% pmask_neg(pmask_neg == 1) = 0;
mpmask_neg(mpmask_neg>=0) = 0;
combined_matrix(:,:,2) = mpmask_neg;

fpmask_pos = fMxM_matrix;
% pmask_pos(pmask_pos == -1) = 0;
fpmask_pos(fpmask_pos<=0) = 0;
combined_matrix(:,:,3) = fpmask_pos;

fpmask_neg = fMxM_matrix;
% pmask_neg(pmask_neg == 1) = 0;
fpmask_neg(fpmask_neg>=0) = 0;
combined_matrix(:,:,4) = fpmask_neg;

%load node definition
path = '/gambit3/fredericks_data/suyeon_data/';
shen = xlsread(sprintf('%s/shen_268_10network_nodecount.xlsx', path));
% shen = xlsread(sprintf('%s/shen_268_labels/shen_268_10network_nodecount.xlsx', path));
x = 3
ind = shen(:,3)==x;
DMNnodesfull = shen(ind,:);
DMNnodes = DMNnodesfull(:,2);

DMN_matrix = combined_matrix(DMNnodes, DMNnodes, :);

save('Facename_DMNedges.mat','DMN_matrix','-v7.3')
