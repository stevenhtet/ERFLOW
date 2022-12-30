EDGESDATA;
NODESDATA;

nodecoords = nodecoords*4/85.667;
nodes = zeros(length(nodecoords),4);
for i = 1:length(nodecoords)
    nodes(i,1) = nodecoords(i,1);
    nodes(i,2) = -1*nodecoords(i,2);
    nodes(i,4) = distance([0 0 0],nodes(i,1:3));
end

edges = zeros(length(edgesdata),3);
for i = 1:length(edgesdata)
    edges(i,1:2) = edgesdata(i,1:2);
    edges(i,3) = distance(nodes(edges(i,1),1:3),nodes(edges(i,2),1:3));
end

edgesinit = edges;
nodesinit = nodes;
E = length(edges);
N = length(nodes);
Einit = length(edgesinit);
Ninit = length(nodesinit);

hold off
figure
for i = 1:length(edgesinit)
    plot([nodes(edgesinit(i,1),1) nodes(edgesinit(i,2),1)],[nodes(edgesinit(i,1),2) nodes(edgesinit(i,2),2)],'.-b')
%     plot3([nodes(edges(i,1),1) nodes(edges(i,2),1)],[nodes(edges(i,1),2) nodes(edges(i,2),2)],[nodes(edges(i,1),3) nodes(edges(i,2),3)],'db')
    hold on
end

% adjmat = genadjmat(edges,N);
adjmatinit = genadjmat(edgesinit,Ninit);
degreevec = zeros(1,Ninit);
for i = 1:Ninit
    degreevec(i) = length(neighbours(i,adjmatinit));
end