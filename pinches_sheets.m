%compute flows in a network due to one contraction-relaxation cycle of a perinuclear sheet;

%sheetT is duration of the contraction (contraction+relaxation lasts
%2*sheetT)
%sheetV is the fluid volume expelled in the contraction




sheetT = 2.5;
sheetV = 10;

Tfactor = 1;
Lfactor = 1;
maxL = 0;
Tend = 1+sheetT*4/4+sheetT*2+2;
h = 0.001*0.1;


energetics = 0;


pinchescell = cell(1,length(1:max(edges(:,2))));
Nsteps = round(Tend/h);
INPUTS = cell(1,Nsteps);
SHEETINPUTS = zeros(1,Nsteps);
QVEC = zeros(E,Nsteps);


SINKS = cell(1,Nsteps); 
    R = 0.03;
    b0 = 0.01*R;
MU = 0.3695846/Tfactor;

M = sinkcount;

%sourcenode L T b bdot

% for v = Ninit+1:Ninit+pinchcount
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
            L = Lmax;
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
        if energetics == 1
            pinchescell{v-Ninit} = [pinchescell{v-Ninit};L T Nstart Nend];
        end

    end

end


%set up sheet contraction here
        starttime = starttimes(i);
        T = Tvec(i);
        Nstart = ceil((1+sheetT/4)/h);
        Nend = ceil((1+sheetT/4+2*sheetT)/h);

        for j = Nstart:Nend
            SHEETINPUTS(j) = sheetV*sin(2*pi*(j-Nstart)/(Nend-Nstart))*pi/2/sheetT;  
        end
        
        if Nend+1<Nsteps
            INPUTS{Nend+1} = [INPUTS{Nend+1}; v, L, T, R, 0, 0];
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


sinkonlynodes = setdiff(sinknodes,sheetnodes);%these are the exit nodes
M1 = length(sinkonlynodes);
M2 = length(sheetnodes);

%compute paths for pressure BCs for exit nodes
sinkpaths = cell(M1-1,1);
sinkpathnum = 0;
for i = 1:M1-1
        sinkpathnum = sinkpathnum + 1;
        sinkpath = findcycle(sinkonlynodes(1),sinkonlynodes(sinkpathnum+1),Tadjmat);
         sinkpaths{sinkpathnum} = sinkpath(1:end-1);
end

%compute paths for pressure BCs for sheet nodes
sheetpaths = cell(M2-1,1);
sinkpathnum = 0;
for i = 1:M2-1
        sinkpathnum = sinkpathnum + 1;
        sinkpath = findcycle(sheetnodes(1),sheetnodes(sinkpathnum+1),Tadjmat);
         sheetpaths{sinkpathnum} = sinkpath(1:end-1);
end



QVEC = zeros(E,Nsteps);
QMAX = zeros(1,Nsteps);

deltan = round((Nend - Nstart)/4);
nstart = Nstart-deltan;
nend = Nend+deltan;
for iii = Nstart:Nend
disp(iii/Nsteps)
inputs = INPUTS{iii};
sheetinput = SHEETINPUTS(iii);
pinchinstance_sheets;%solve system instantaneously
QVEC(:,iii) = qvec;%these are the computed fluxes
QMAX(iii) = qmax;
end

