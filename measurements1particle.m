edgeltv = zeros(length(nodetime)/2,2);
for j = 1:length(nodetime)/2
        v = nodetime(2*j-1,1); w = nodetime(2*j,1);
        t = (nodetime(2*j,2) - nodetime(2*j-1,2))*h;
        e = edgesinit2scalar(min(v,w),max(v,w));
        l = edgesinit(e,3);
        edgeltv(j,:) = [e l/t];
end

numvisited = length(nodetime/2);
if ~isempty(nodetime)
numvisitednorep = length(setdiff(nodetime(:,1),[]));
else
    numvisitednorep = 1;
end
