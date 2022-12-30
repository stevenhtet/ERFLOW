if isempty(inputs) %if no pinches 
    qvec = zeros(E,1);
    qmax = 0;
else
    sourcenodes = inputs(:,1);
    Lvec = inputs(:,2);
    Tvec = inputs(:,3);
    bvec = inputs(:,4);
    bdotvec = inputs(:,5);
    
    %modify pressure drop matrix with Hagen-Poiseuille for pinching tubules
    %(pinch nodes previously termed source nodes)
    for k = 1:length(sourcenodes)
        node = sourcenodes(k);
        L = Lvec(k);
        b = bvec(k);
        list = neighbours(node,adjmat);
        for v = list
            if b~=R
            l1 = distance(nodes(node,1:3),nodes(v,1:3)) - L;
            chimat(node, v) = 8*mu/pi*(l1/R^4 + (1/b^3 - 1/R^3)*L/3/(R-b));
            chimat(v, node) = chimat(node,v);
            else
            chimat(node,v) = gamma*distance(nodes(node,1:3),nodes(v,1:3));
            chimat(v,node) = chimat(node,v);
            end
        end
    end



    %set up K2 eqns
    PMAT = zeros(E-N+1,E+M);
    for i = 1:E-N+1
        cycle = cycles{i};
        vec = zeros(1,E);
        for j = 1:length(cycle)-1
            v = cycle(j);
            w = cycle(j+1);
            if v<w
                vec(edges2qmat(v,w)) = vec(edges2qmat(v,w))+1*chimat(v,w);
            else
                vec(edges2qmat(w,v)) = vec(edges2qmat(w,v))-1*chimat(w,v);                
            end
        end
        PMAT(i,1:E) = vec;
    end
    
    %K1 at exit nodes; 'sinks' suffix refers to exit nodes,
    %previously termed sink nodes
    MMATsinks = zeros(N,M);
    for i = 1:M
        MMATsinks(sinknodes(i),i) = -1;
    end
    MASSMAT = [MMAT MMATsinks];
    
   
    %set up pressure BC equations; 
    PMATsinks = zeros(M-1,E+M);
    for i = 1:M-1
        cycle = sinkpaths{i};
        vec = zeros(1,E);
        for j = 1:length(cycle)-1
            v = cycle(j);
            w = cycle(j+1);
            if v<w
                vec(edges2qmat(v,w)) = vec(edges2qmat(v,w))+1*chimat(v,w);               
            else
                vec(edges2qmat(w,v)) = vec(edges2qmat(w,v))-1*chimat(w,v);               
            end
        end
        PMATsinks(i,1:E) = vec;
    end
    
    
    
    %put together K1, K2 and pressure BC eqns into one matrix
    Qvec = inputs(:,6);
    MAT = [MASSMAT; PMAT; PMATsinks];
    VEC = zeros(E+M,1);
    VEC(sourcenodes) = -Qvec;

    soln = mldivide(MAT,-VEC); %invert linear system for soln = [fluxes through edges; sources at exit nodes]'
    qvec = soln(1:E);
    qmax = max(abs(qvec));
    vmax = qmax/pi/R^2;
    if sum(soln(E+1:end)) + sum(Qvec) > 10^-10
        disp('error')
    end
end
         