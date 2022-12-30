function Tedges = gentree(dwmatadj)
%BFS algorithm to generate spanning tree from adjacency matrix (distance-weighted or otherwise)
Tnodes = [1];
Tedges = [];
unexp = [1];
while ~isempty(unexp)
    v = unexp(1);
    unexp(1) = [];
    list = find(dwmatadj(v,:)~=0);
    for i = 1:length(list)
        w = list(i);
        if isempty(find(unexp == w, 1)) && isempty(find(Tnodes == w, 1))%(v,w) doesn't produce cycle
            unexp = [unexp w];
            Tedges = [Tedges; v w dwmatadj(v,w)];
            Tnodes = [Tnodes w];
        end
    end
end