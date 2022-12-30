function CYCLE = findcycle(v0,v1,Tadjmat)
%finds a path from v0 to v1 on a spanning tree with adjcency matrix Tadjmat
yess = 0;
explored = [];
queue = cell(1);
queue{1} = [v0];
while ~isempty(queue)
    path = queue{1};
    queue(1) = [];
    node = path(end);
    if isempty(find(explored==node,1))
        neighbors = find(Tadjmat(node,:)==1);
        for j = 1:length(neighbors)
            if isempty(find(explored==neighbors(j),1))
            newpath = [path neighbors(j)];
            queue{end+1} = newpath;
            if neighbors(j) == v1
                CYCLE = newpath;
                queue = [];
            end
            end
        end
        explored = [explored node];   
    end
end
CYCLE = [CYCLE v0];
    
            
%      4   233   104   381   174   493   174     4
    


