%construct honeycomb network

LL = 1; %set edge length
up = [sqrt(3)/2 1/2 0 0]*LL;
down = [sqrt(3)/2 -1/2 0 0]*LL;
nodes = [];
Nhex = 3;   %size of network

for i = 1:Nhex
    nodes = [nodes; [0 LL 0 0]*(i+(ceil(i/2)-1))];
   if i == Nhex
       nodes(end,4) = 1;
   end
   for j = 1:2*Nhex-i 
       if mod(i+j,2)==0
            nodes = [nodes; nodes(end,:)+down];
       else
           nodes = [nodes; nodes(end,:)+up];
       end
       if j == 2*Nhex-i || j == 2*Nhex-i-1 
           nodes(end,4) = 1;
       end
   end
end

nn = length(nodes);
for i = 1:nn
    if nodes(i,1) ~= 0
        nodes = [nodes; nodes(i,1)*-1 nodes(i,2) nodes(i,3) nodes(i,4)];
    end
end

nodes = [nodes; nodes(:,1) nodes(:,2)*-1 nodes(:,3) nodes(:,4)];


sinknodes = [];
for i = 1:length(nodes)
    if nodes(i,4) ~= 0  
        sinknodes = [sinknodes; i];
    end
end
Ninit = length(nodes);
N = Ninit;

Dmat = zeros(N,N);
for i = 1:N
    for j = i+1:N
        Dmat(i,j) = distance(nodes(i,1:3),nodes(j,1:3));
    end
end
Dmat = Dmat + Dmat';
for i = 1:N
    Dmat(i,i) = NaN;
end
edges = [];

for i = 1:N
    for j = i+1:N
        if Dmat(i,j) < 1.01*LL
            edges = [edges; i j LL];
        end
    end
end



Einit = length(edges);
E = Einit;
adjmatinit = genadjmat(edges,N);
    pinchcount = E;
edgesinit = edges;
nodesinit = nodes;
edgelengths = edgesinit(:,3);
pinchvec = 1:E;
hold off
pinchnodes = zeros(pinchcount,4);
pinchedges = zeros(pinchcount*2,3);
for i = 1:length(pinchvec)
    pinchnodes(i,:) = (nodes(edges(pinchvec(i),1),:)+nodes(edges(pinchvec(i),2),:))/2; %can maake more general
    pinchnodes(i,4) = distance(pinchnodes(i,1:3),[0 0 0]);
    plot(pinchnodes(i,1),pinchnodes(i,2),'xr')
    hold on
    pinchedges(2*i-1,:) = [edges(pinchvec(i),1) N+i distance(nodes(edges(pinchvec(i),1),1:3),pinchnodes(i,1:3))];    
    pinchedges(2*i,:) = [edges(pinchvec(i),2) N+i distance(nodes(edges(pinchvec(i),2),1:3),pinchnodes(i,1:3))];
end
    edges(pinchvec,:)=[];
nodes = [nodes; pinchnodes];
edges = [edges; pinchedges];
N = length(nodes);
Nsources = N;
Esources = length(edges);

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

%plot graphs
for i = 1:length(Tedges)
    plot([nodes(Tedges(i,1),1) nodes(Tedges(i,2),1)],[nodes(Tedges(i,1),2) nodes(Tedges(i,2),2)],'.-k')
end
for i = 1:length(edges)-1
    plot([nodes(edges(i,1),1) nodes(edges(i,2),1)],[nodes(edges(i,1),2) nodes(edges(i,2),2)],'.-b');hold on
end

sinkcount = length(sinknodes);
randsink = 0;





    