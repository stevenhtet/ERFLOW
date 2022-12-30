edgevels = zeros(1,Einit);
edgenums = zeros(1,Einit);
ES = EV(:,1);
VS = EV(:,2);
for i = 1:Einit
    edgevels(i) = mean(VS(ES == i));
    edgenums(i) = length(find(ES==i));
end



MSD = mean(DISPS,2);

if NOFLOW == 0
edgevelsNOFLOW = zeros(1,Einit);
edgenumsNOFLOW = zeros(1,Einit);
ESNOFLOW = EVNOFLOW(:,1);
VSNOFLOW = EVNOFLOW(:,2);
for i = 1:Einit
    edgevelsNOFLOW(i) = mean(VSNOFLOW(ESNOFLOW == i));
    edgenumsNOFLOW(i) = length(find(ESNOFLOW==i));
end
MSDNOFLOW = mean(DISPSNOFLOW,2);
end

TS = ES;
for i = 1:length(EV)
  TS(i) = edges(ES(i),3)/VS(i);
end
    
    
    