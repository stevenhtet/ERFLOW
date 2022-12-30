%compute flows in a network with pinching tubules;
%tubules pinch Tfactor times faster and more frequently than Holcman 2018
%parameters
%If Lmax = 1, pinch length = (L - 2R)*LFACTOR
%If lmax = 1, pinches are drawn from a normal distribution and are Lfactor
%times the pinches in Holcman 2018

Tfactor = 10;
Lfactor = 1;
LFACTOR = 1;
maxL = 1;
Tend = 6;
h = 0.001*0.05;

pinchescell = cell(1,length(Ninit+1:max(edges(:,2))));
Nsteps = round(Tend/h);
INPUTS = cell(1,Nsteps);
QVEC = zeros(E,Nsteps);


SINKS = cell(1,Nsteps); 
    R = 0.03;
    b0 = 0.01*R;
MU = 0.3695846/Tfactor;


M = sinkcount;

%set up pinches on each edge; here v loops through nodes with index > Ninit
%i.e. the pinch nodes
for v = Ninit+1:max(edges(:,2))
    starttimes = [exprnd(MU)];
    Tvec = [getT(Tfactor)];
    while starttimes(end) < Tend
        Tint = exprnd(MU);
        starttimes = [starttimes, starttimes(end) + 2*Tvec(end) + Tint];
        Tvec = [Tvec, getT(Tfactor)];
    end
    Lmax = (sum(edges(find(edges(:,2)==v),3)))/2 - 2*R;
        
    for i = 1:length(starttimes) - 2
        starttime = starttimes(i);
        if maxL == 1
            L = Lmax*LFACTOR;
        else
            L = getL()*Lfactor;
            if L>Lmax
                L = Lmax;
            end
        end
        T = Tvec(i);
        Nstart = ceil(starttime/h);
        Nend = Nstart + ceil(2*T/h);

        check = 0;
        for j = Nstart:min(Nend,Nsteps)
            b = (R+b0)/2 + (R-b0)/2*cos(2*pi*(j-Nstart)/(Nend-Nstart));
            bdot = -(R-b0)*pi/2/T*sin(2*pi*(j-Nstart)/(Nend-Nstart));
           Qj = -2*pi*L*bdot*(R - 2*(R-b)/3);
            INPUTS{j} = [INPUTS{j}; v, L, T, b, bdot, Qj];  
        end
        
        if Nend+1<Nsteps
        INPUTS{Nend+1} = [INPUTS{Nend+1}; v, L, T, R, 0, 0];
        end
    end          
end


%set up matrix of pressure drop coeffs chi across each edge
chimat = zeros(N,N);
mu = 10^-9;
gamma = 8*mu/pi/R^4;
for i = 1:E
    chimat(edges(i,1),edges(i,2)) =  gamma*edges(i,3);    
end
chimat = chimat + chimat';
chimatreset = chimat;

%set up K1 equations
MMAT = zeros(N,E);
for v = 1:N
    list = neighbours(v,adjmat);
    vec = zeros(1,E);
    for i = 1:length(list)
        w = list(i);
        if v<w
            vec(edges2qmat(v,w)) = 1;
        else
            vec(edges2qmat(w,v)) = -1;
        end
    end
    MMAT(v,:) = vec;
end


%generate paths between exit node 1 and exit node 2 to M
%(findcycle.m in hindsight is more aptly named findpath.m...)
sinkpaths = cell(M-1,1);
sinkpathnum = 0;
for i = 1:M-1
        sinkpathnum = sinkpathnum + 1;
        sinkpath = findcycle(sinknodes(1),sinknodes(sinkpathnum+1),Tadjmat);
         sinkpaths{sinkpathnum} = sinkpath(1:end-1);
end



QVEC = zeros(E,Nsteps);
QMAX = zeros(1,Nsteps);
nstart = round(1 + 1/Tend*Nsteps);
nend = round(Nsteps - 1/Tend*Nsteps);
for iii = nstart:nend
disp(iii/Nsteps)
inputs = INPUTS{iii};
pinchinstance_pressureBC; %solve system instantaneously
QVEC(:,iii) = qvec; %these are the computed fluxes
QMAX(iii) = qmax;
end

