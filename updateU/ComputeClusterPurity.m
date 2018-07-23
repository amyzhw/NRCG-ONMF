function [purity, cluster_purities, T, std_purity] = ComputeClusterPurity(IDX, C, k, nc)
if ~exist('nc', 'var') || isempty(nc)
    nc = k;
end
n = size(IDX,1);

classes = clabel2dataclasses(C,nc);
clusters = clabel2dataclasses(IDX,k);

T = clusters' * classes;

purity = full(sum(max(T, [], 2)) / n);

if nargout>1
    cluster_purities = zeros(nc,1);
    indnz = sum(T,2) > 0;
    cluster_purities(indnz) = max(T(indnz,:),[],2)./sum(T(indnz,:),2);
    weights = sum(T,2) / n;
    std_purity = sqrt(weights'*(cluster_purities - mean(cluster_purities)).^2);  
end