function ratio = association( conn_mat, sigma, label)
% conn_mat is N by N matrix
% label is 1D array of network labels

no_net = max(label);
no_nodes = size( conn_mat, 1);
ratio = zeros(1, no_net);

% convert the correlation matrix to weight matrix
dist = (1-conn_mat);
weight = exp(-dist/sigma);
% set the diagonal of the weight matrix to be 0
for i = 1:no_nodes
    weight(i, i) = 0;
end

for nn = 1:no_net;
    cur_id = find( label==nn);
    link_within = sum(sum( weight(cur_id, cur_id)));
    link_across = sum( sum( weight(cur_id, :)));
    ratio(nn) = link_within/link_across;
end
return;