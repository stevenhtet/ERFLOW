figure
subplot(1,3,1)

AAA = INSTVELS(2:end-2,:);
hold on; histogram(AAA(AAA~=0),'Normalization','pdf','DisplayStyle','stairs','EdgeColor','b')
xlabel('instantaneous speeds (μm/s)');
subplot(1,3,2)
hold on; histogram(VS,'Normalization','pdf','DisplayStyle','stairs','EdgeColor','b')
xlabel('edge traversal speeds (μm/s)');
subplot(1,3,3)
hold on; histogram(edgevels ,'Normalization','pdf','DisplayStyle','stairs','EdgeColor','b')
xlabel('average edge traversal speeds (μm/s)');
