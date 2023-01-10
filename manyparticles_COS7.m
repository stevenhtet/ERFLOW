%track Brownian particles in network subject to computed flows
%total #particles = nump*#nodes in original graph

%set NOFLOW = 1 if pure diffusion; still need to run part of
%pinches_pressureBC.m to define INPUTS cell and UVEC array

nump = 1;
D = 0.6; %diffusivity
sigmar = (2*D*h)^0.5; 
rparticle = 0.0025;%
trecord = 0.018; %time resolution for 'instantaneous' speeds
nrecord = round(trecord/h);
maxnodetime = h;
NOFLOW = 0;

%initialise
hex = 0;
Reff = R - rparticle;
rval = @(x)(x(1)^2 + x(2)^2)^0.5;
nstart = 1 + 1/Tend*Nsteps;
nend = Nsteps - 1/Tend*Nsteps;
UVEC = QVEC(:,nstart:nend)/pi/R^2;
clear('QVEC');
INPUTS(nend+1:end)=[];
INPUTS(1:nstart-1)=[];
edgesinit(:,1:2) = sort(edgesinit(:,1:2),2);
edgesinit2scalar = zeros(Ninit,Ninit);
for i = 1:Einit
    edgesinit2scalar(edgesinit(i,1),edgesinit(i,2)) = i;
    edgesinit2scalar(edgesinit(i,2),edgesinit(i,1)) = i;
end
if edges(end,1)~=0
edges = [edges; 0 0 0];
end

zvectors = zeros(E,2);
xvectors = zeros(E,2);
for i = 1:E
    A = nodes(edges(i,1),1:2);
    B = nodes(edges(i,2),1:2);
    zhat = B - A;
    zhat = zhat/(zhat(1)^2+zhat(2)^2)^0.5;
%     xhat = -y x
    xhat = [-zhat(2) zhat(1)];
    if A(1)>B(1)
        xhat = -xhat;
    elseif A(1)==B(1) && A(2)>B(2)        
        xhat = -xhat;
    end
    xvectors(i,:) = xhat;
    zvectors(i,:) = zhat;
end

    FLOW = 1;%this is always 1; ignore confusion with variable NOFLOW
    meandisps = zeros(round(length(UVEC)/nrecord),1);
    meannumvisited = 0;
    meannumvisitednorep = 0;
    NUMV = [];
    NUMVnorep = [];
    EV = [];

    nonsinknodes = setdiff(1:Ninit,sinknodes);
    numparticles = nump*length(nonsinknodes);
    DISPS = zeros(round(length(UVEC)/nrecord),numparticles);
    INSTVELS = zeros(round(length(UVEC)/nrecord),numparticles);
    
    %if pure diffusion remove INPUTS
    if NOFLOW == 1
        for inputindex = 1:length(INPUTS)
            INPUTS{inputindex} = [];
        end
        UVEC = UVEC*0;
    end
    %main loop
    for aaaa = 1:length(nonsinknodes)
        for repeat = 1:nump
            startnode = nonsinknodes(aaaa);
            aaaa
            track1particle;
            measurements1particle; %compute edge traversal speeds
            DISPS(:,aaaa) = disps;
            INSTVELS(:,aaaa) = velocities;
            EV = [EV; edgeltv]; %list of [edge number, edge traversal speed] for every edge traversal event
            NUMV = [NUMV; numvisited];
            NUMVnorep = [NUMVnorep; numvisitednorep];   
        end
    end
  
    
    %calculate average edge traversal speeds
    edgevels = zeros(1,Einit);
    edgenums = zeros(1,Einit);
    ES = EV(:,1);
    VS = EV(:,2);
    for i = 1:Einit
        edgevels(i) = mean(VS(ES == i));
        edgenums(i) = length(find(ES==i));
    end

    meaninstvels = mean(INSTVELS(:));

    
