%%Construct C1-4 networks.

C1DATA;
%%specify that we want to construct C1 network; use C2DATA, C3DATA, C4DATA
%%instead for C2-C4


twopinchthreshold = 100;
pinchthreshold = 0.1;

%construct pinch point
edgelengths = edgesinit(:,3);
twopinchvec = find(edgelengths>twopinchthreshold);
midpinchvec = [];
for i = 1:length(twopinchvec)
    p = rand;
    if p>0.5
        midpinchvec = [midpinchvec; twopinchvec(i)];
    end
end
hold off;
midpinchcount = length(midpinchvec);
midpinchvec = sort(midpinchvec);
midpinchnodes = zeros(midpinchcount,4);
midpinchedges = zeros(midpinchcount*2,3);
for i = 1:length(midpinchvec)
    elength = distance(nodes(edges(midpinchvec(i),1),1:3),nodes(edges(midpinchvec(i),2),1:3));
    alpha = 0.5;
    midpinchnodes(i,:) = alpha*(nodes(edges(midpinchvec(i),1),:))+(1-alpha)*nodes(edges(midpinchvec(i),2),:); %can maake more general
    midpinchnodes(i,4) = distance(midpinchnodes(i,1:3),[0 0 0]);

    midpinchedges(2*i-1,:) = [edges(midpinchvec(i),1) N+i distance(nodes(edges(midpinchvec(i),1),1:3),midpinchnodes(i,1:3))];    
    midpinchedges(2*i,:) = [edges(midpinchvec(i),2) N+i distance(nodes(edges(midpinchvec(i),2),1:3),midpinchnodes(i,1:3))];
end
    edges(midpinchvec,:)=[];
nodes = [nodes; midpinchnodes];
edges = [edges; midpinchedges];
N = length(nodes);
Nsources = N;
Esources = length(edges);


%construct pinch point
edgelengths = edges(1:length(edges)-length(midpinchedges),3);
pinchvec = find(edgelengths>pinchthreshold);
pinchvec = [pinchvec;(length(edges)-length(midpinchedges)+1:length(edges))'];
pinchcount = length(pinchvec);
pinchvec = sort(pinchvec);
pinchnodes = zeros(pinchcount,4);
pinchedges = zeros(pinchcount*2,3);
for i = 1:length(pinchvec)
    elength = distance(nodes(edges(pinchvec(i),1),1:3),nodes(edges(pinchvec(i),2),1:3));
    alpha = 0.5;
    pinchnodes(i,:) = alpha*nodes(edges(pinchvec(i),1),:)+(1-alpha)*nodes(edges(pinchvec(i),2),:); %can maake more general
    pinchnodes(i,4) = distance(pinchnodes(i,1:3),[0 0 0]);
    pinchedges(2*i-1,:) = [edges(pinchvec(i),1) N+i distance(nodes(edges(pinchvec(i),1),1:3),pinchnodes(i,1:3))];    
    pinchedges(2*i,:) = [edges(pinchvec(i),2) N+i distance(nodes(edges(pinchvec(i),2),1:3),pinchnodes(i,1:3))];
end
    edges(pinchvec,:)=[];
nodes = [nodes; pinchnodes];
edges = [edges; pinchedges];
N = length(nodes);
Nsources = N;
Esources = length(edges);


sinkcount = length(sinknodes);

N = length(nodes);
E = length(edges);


% construct adjacency matrix
dwadjmat = gendwadjmat(edges,N);
adjmat = genadjmat(edges,N);

%construct spanning tree
Tedges = gentree(dwadjmat);



%generate cycle basis
edges(:,1:2) = sort(edges(:,1:2),2);
Tedges(:,1:2) = sort(Tedges(:,1:2),2);
Tadjmat = genadjmat(Tedges,N);
cycles = cell(E-N+1,1);
cyclenum = 0;
basisedges = setdiff(edges,Tedges,'rows');
for i = 1:E-N+1
        cyclenum = cyclenum + 1;
        cycles{cyclenum} = findcycle(basisedges(i,1),basisedges(i,2),Tadjmat);
end





%enumerate/define flow rate variables of edges
E = length(edges);
edges2qmat = zeros(N,N);
for i = 1:E
    edges2qmat(edges(i,1),edges(i,2)) = i;
end

R = 0.03;
b0 = 0.01*R;


for i = 1:length(edgesinit)
    hold on
    plot([nodesinit(edgesinit(i,1),1) nodesinit(edgesinit(i,2),1)],[nodesinit(edgesinit(i,1),2) nodesinit(edgesinit(i,2),2)],'.-k')
end

