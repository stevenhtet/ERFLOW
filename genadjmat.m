function mat = genadjmat(edges,N)
%generate adjancency matrix from list of edges; 
%N = number of edges;
%edges = edges (i,j) listed as an Nx2 matrix
mat = zeros(N,N);
for i = 1:length(edges)
    mat(edges(i,1),edges(i,2)) = 1;
end
mat = mat+ mat';