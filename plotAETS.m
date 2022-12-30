figure(1)
H = histogram(edgevels,'Normalization','pdf');
        binedges = H.BinEdges;
        for i = 1:(length(binedges)-1)
            xx(i) = binedges(i) + binedges(i+1);
        end
        xx = xx/2;
        yy = H.Values;

figure(2)
subplot(1,2,plotindex)
plot(xx(1:length(yy)),yy,'ob','MarkerFaceColor','b');

xlabel('average edge traversal speeds (μm/s)');

DIST = fitdist(edgevels','normal');
range = binedges(end); 
x = 0:range/250:range;
y = pdf('normal',x,DIST.mu,DIST.sigma);
hold on;
plot(x,y,['b-'])
legend(['mean = ' num2str(mean(edgevels(~isnan(edgevels)))) 'μm/s'],['mean = ' num2str(DIST.mu) 'μm/s']) 