%simulates the motion of one brownnian particle in network subject to
%computed flows; particle reflects of tubule walls

nodetime = []; %[node index, time] concatanated to this list at the beginning and end of each edge traversal;
P = zeros(length(UVEC),5); %each row is the position vector of the particle 
%at a given time in the form [edge, node, x, y, z]; edge = 0 when particle
%is in node, node = 0 when particle is in edge;
%z is coordinate along tubule, z = 0 at beginning of edge and z = edge
%length at end of edge
%for an edge (i,j), we consider the smaller of i and j to be the 'beginning' of the edge.

%initialise particle in given node
theta = 2*pi*rand;
R0 = Reff*sqrt(rand);
P(1,:) = [0 startnode R0*cos(theta) R0*sin(theta) edges(1,3)*rand];

p = rand;
pass = 1;
prob = 0;
disps = zeros(round(length(UVEC)/nrecord),1);
positions = zeros(round(length(UVEC)/nrecord),3);
velocities = zeros(round(length(UVEC)/nrecord),1);
modcount = 2;
inpinch = 0;
b = R;
dvec = zeros((length(UVEC))-1,3);
bvec = dvec;
outside = 0;
for n = 1:(length(UVEC))-1
    %record particle position
    if mod(n,nrecord)==0 && modcount<round(length(UVEC)/nrecord)
        if P(n,2) ~= 0 %particle in node
            disps(modcount) = distance(nodes(startnode,1:3),nodes(P(n,2),1:3))^2;
            positions(modcount,:) = nodes(P(n,2),1:3);           
        else
           eee = P(n,1);
           xhat = xvectors(eee,:);
           vvv = edges(eee,1);
           www = edges(eee,2);
           lll = edges(eee,3);
           zzz = P(n,5);
           aaa = nodes(vvv,1:3); bbb = nodes(www,1:3);
           ccc = aaa + zzz/lll*(bbb-aaa);
           disps(modcount) = distance(nodes(startnode,1:3),ccc)^2;
            positions(modcount,:) = ccc+P(n,3)*[xhat 0];
        end
        modcount = modcount+1;
    end
    
    e = P(n,1);
    if e ~= 0 %if particle in edge
       %CHECK IF PARTICLE IN PINCH
        potentialpinchnode = max(edges(e,2));
        if ~isempty(INPUTS{n})  %INPUTS is cell of currently active pinches
            pinchnodeindex = find(INPUTS{n}(:,1)==potentialpinchnode,1);
        else
            pinchnodeindex = [];
        end
        
        if isempty(pinchnodeindex)
            inpinch = 0;
            pinchpresent = 0;
        else
            pinchpresent = 1;
            L = INPUTS{n}(pinchnodeindex,2);
            b = INPUTS{n}(pinchnodeindex,4);
            bdot = INPUTS{n}(pinchnodeindex,5);
        end
        if edges(e,3)-P(n,5)>L
            inpinch = 0;
        else
            inpinch = 1;
        end
        
        %x0 is spatial location of particle
        %if particle in pinch, and if pinch has advanced since previous timestep
        %such that particle now lies outside of pinch, artificially
        %move particle back inside pinch
        x0 = P(n,3:5); 
        r = rval(P(n,3:5));
        if pinchpresent==0
            inpinch = 0;
        else
            if edges(e,3)-x0(3)>L %if x0 not in pinch
                inpinch = 0;
            else 
                inpinch = 1;
                if R~=b
                    xi = L*b/(R-b);                  
                    a0 = b + b/xi*(edges(e,3)-x0(3))-rparticle; 
                    if r>a0
                        P(n,3:4) = P(n,3:4)*a0/r;
                        x0 = P(n,3:5);
                        r = rval(x0);
                    end
                else
                    if r>Reff
                    0;
                    end
                end
            end
        end
        
        
            %update particle position according to advection and diffusion
            if inpinch == 0
                P(n+1,:) = P(n,:) + [0 0 0 0 UVEC(e,n)*2*(1-r^2/R^2)*h] + [0 0 normrnd(0,sigmar) normrnd(0,sigmar) normrnd(0,sigmar)];
            else
                dadt = (P(n,5)-(edges(e,3)-L))/L*bdot;
                dadz = -(R-b)/L;
                a = R - (R-b)*(P(n,5)-(edges(e,3)-L))/L;
                x = P(n,3);
                y = P(n,4);
                r = rval(P(n,3:5));
                U = (UVEC(e,n)*pi*R^2 - 2*pi*bdot/L*(P(n,5)-(edges(e,3)-L))^2*(R/2-(R-b)*(P(n,5)-(edges(e,3)-L))/3/L))/pi/a^2;
                V = dadt*r/a*(2-(r/a)^2)+2*dadz*U*(r/a)*(1-(r/a)^2);
                
                P(n+1,:) = P(n,:) + [0 0 x/r*V*h y/r*V*h U*2*(1-r^2/R^2)*h] + [0 0 normrnd(0,sigmar) normrnd(0,sigmar) normrnd(0,sigmar)];

            end


        %
        d = P(n+1,3:5)-P(n,3:5);
        xfinal = P(n+1,3:5);
        r = sqrt(xfinal(1)^2+xfinal(2)^2);
        %calculate tubule radius aa at longitudinatl coord or particle according to whether particle is in pinch 
        if pinchpresent==0
            aa = Reff;
            inpinch = 0;
        else
            if edges(e,3)-xfinal(3)>L %xfinal not in pinch
            aa = Reff;
            inpinch = 0;
            else %xfinal in pinch
                if R~=b
                    xi = L*b/(R-b);                  
                    aa = b + b/xi*(edges(e,3)-xfinal(3))-rparticle; 
                    inpinch = 1;
                else
                    inpinch = 0;
                    aa = Reff;
                end
            end
        end
        whilecounter = 0;
        %reflect at wall if updated particle position is outside tubules
        while r>aa && whilecounter<1000
            updatex = 1;
            if inpinch == 0 %reflect at cylindrical wall
                A = d(1)^2 + d(2)^2; 
                B = 2*(d(1)*x0(1)+d(2)*x0(2));
                C = x0(1)^2+x0(2)^2-(Reff)^2;
                alpha = (-B+(B^2-4*A*C)^0.5)/2/A;

                xp = x0 + alpha*d;
                rp = (xp(1)^2+xp(2)^2)^0.5;
                nvec = [xp(1)/rp xp(2)/rp 0];
                xfinal = x0 + d - 2*(1-alpha)*dot(nvec,d)*nvec;
                r = sqrt(xfinal(1)^2+xfinal(2)^2);
                aa = Reff;
                x0bk = x0;
                xpbk = xp;
                xfinalbk = xfinal;
                dbk = d;               
                if pinchpresent == 1 && b~=R
                    if edges(e,3)-xp(3)<L
                        inpinch = 1;
                        updatex = 0;                       
                    end
                end     
            else %inpinch == 1, reflect at pinching wall
                if P(n+1,1)~=P(n,1)
                    x0(3) = edges(e,3)+(edges(P(n-1,1),3)-x0(3));
                end                        
                        le = edges(e,3);
                        if R~=b
                        xi = L*b/(R-b);
                        else
                        inpinch = 0;
                        updatex = 0;
                        whilecounter = whilecounter+1;
                        end
                        zeta = b - rparticle + b*le/xi-b*x0(3)/xi;
                        A = d(1)^2 + d(2)^2 - d(3)^2*b^2/xi^2;
                        B = 2*(d(1)*x0(1)+d(2)*x0(2)+d(3)*zeta*b/xi);
                        C = x0(1)^2+x0(2)^2-zeta^2;
                        alpha = (-B+(B^2-4*A*C)^0.5)/2/A;

                        costheta = (xi+L)/((xi+L)^2+R^2)^0.5;
                        sintheta = R/((xi+L)^2+R^2)^0.5;
                        xp = x0 + alpha*d;

                        rp = (xp(1)^2+xp(2)^2)^0.5;
                        nvec = [costheta*xp(1)/rp costheta*xp(2)/rp sintheta];
                        xfinal = x0 + d - 2*(1-alpha)*dot(nvec,d)*nvec;                

                %calculate r and aa
                    r = sqrt(xfinal(1)^2+xfinal(2)^2);
                    if R~=b
                    if xfinal(3)>le
                        aa = b - b/xi*(edges(e,3)-xfinal(3))-rparticle;
                    elseif le - xfinal(3)>L
                        aa = R - rparticle;
                        inpinch = 0;
                    else                    
                    aa = b + b/xi*(edges(e,3)-xfinal(3))-rparticle; 
                    end
                    else
                    aa = Reff;
                    end
                    if le - xp(3)>L
                        inpinch = 0;
                        updatex = 0;
                        whilecounter = whilecounter+1;
                    end    
                    x0bk = x0;
                    xfinalbk = xfinal;
                    xpbk = xp;
                    dbk = d;
            end
            if updatex == 1
            x0 = xp;
            d = xfinal - x0;
            whilecounter = whilecounter+1;
            end
        end
        if whilecounter>999
            pass = 0; 
        end
        lastL = L;
        lastb = b;
        P(n+1,3:5) = xfinal;
        
        
        %reach end of edge
        if P(n+1,5)>edges(e,3)
            if edges(e,2)>Ninit  %pass through pinchnode
                list = neighbours(edges(e,2),adjmat);
                vnew = list(list~=edges(e,1));
                enew = edges2qmat(min(edges(e,2),vnew),max(edges(e,2),vnew));
%                 pass = 0;
                P(n+1,1) = enew;
                P(n+1,5) = edges(enew,3) - (P(n+1,5) - edges(e,3));
            else %normal node
                if edges(e,2) ~=  lastnodetime(1)
                    nodetime = [nodetime; lastnodetime; edges(e,2) n+1];
                end
                P(n+1,:) = [0 edges(e,2) 0 0 0];
                p = rand;
            end
        end
        
        %reach beginning of edge; 
        if P(n+1,5)<0
            if edges(e,1)>Ninit %pass through pinch node
                list = neighbours(edges(e,1),adjmat);
                vnew = list(list~=edges(e,2));
                enew = edges2qmat(min(edges(e,1),vnew),max(edges(e,1),vnew));
                P(n+1,1) = enew;
                P(n+1,5) = -P(n+1,5);
            else
                if edges(e,1)~= lastnodetime(1)
                nodetime = [nodetime; lastnodetime; edges(e,1) n+1];
                end
                P(n+1,:) = [0 edges(e,1) 0 0 0];
                p = rand;
            end
        end
    elseif P(n,2)~= 0 %particle is in node
        if p<=prob %escape node if p<prob
            if FLOW == 1
                for k = 1:length(listt)
                    w = listt(k);
                    if w<v
                        probs(k) = probs(k) - UVEC(edges2qmat(w,v),n); %#ok<SAGROW>
                    else
                        probs(k) = probs(k) + UVEC(edges2qmat(v,w),n); %#ok<SAGROW>
                    end
                end
                for k = 1:length(probs)
                    probs(k) = 1 + probs(k)*R/D;
                end
                probs(probs<=0) = 0;
                psum = sum(probs);
                if psum~=0
                    probs = probs/psum;
                    probs = cumsum(probs);
                    probss = [0 probs'];
                    pp = rand;
                    for k = 1:length(probs)
                        if pp<=probss(k+1) && pp>probss(k)
                            vv = listt(k);
                        end
                    end
                else
                    vv = listt(randi([1 length(listt)]));
                end
            else
                vv = listt(randi([1 length(listt)]));
            end

            %move particle into edge
            if v<vv
                e = edges2qmat(v,vv);
                if n>2
               P(n+1,:) = [e, 0, P(n-2,3), P(n-2,4), 0];
                else
               P(n+1,:) = [e, 0, P(n-1,3), P(n-1,4), 0];
                end
               prob = 0;
            else
                e = edges2qmat(vv,v);
               if n>2
                 P(n+1,:) = [e, 0, P(n-2,3), P(n-2,4), edges(e,3)];
                else
                    P(n+1,:) = [e, 0, P(n-1,3), P(n-1,4), edges(e,3)];
                end
               prob = 0;
            end
            lastnodetime = [v n+1];%enter edge
        else
            if prob == 0
            v = P(n,2);
            listt = neighbours(v,adjmat);
            probs = zeros(length(listt),1);                
            end
            if FLOW == 1
                for k = 1:length(listt)
                    w = listt(k);
                    if w<v
                        probs(k) = probs(k) - UVEC(edges2qmat(w,v),n);
                    else
                        probs(k) = probs(k) + UVEC(edges2qmat(v,w),n);
                    end
                end
                for k = 1:length(probs)
                    probs(k) = 1 + probs(k)*R/D;
                end
            end
            prob = prob + 1/(maxnodetime/h);            
            P(n+1,2) = P(n,2);
        end
    end

end

%'instantaneous' speeds
for i = 1:length(positions)-2
    velocities(i) = distance(positions(i+1,:),positions(i,:))/trecord; 
end

rvec = (P(:,3).^2+P(:,4).^2).^0.5;
evec = P(:,1);
