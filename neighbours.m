function list = neighbours(v,adjmat)
list = find(adjmat(v,:)==1);