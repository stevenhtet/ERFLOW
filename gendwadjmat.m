function mat = gendwadjmat(edges,N)
%generate distance-weighted adjancency matrix from list of edges; 
%N = number of edges;
%edges = edges (i,j) listed as an Nx2 matrix
mat = zeros(N,N);
for i = 1:length(edges)
    mat(edges(i,1),edges(i,2)) = edges(i,3);
end
mat = mat + mat';