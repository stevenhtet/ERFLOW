% function e = findedge(x1,y1,x2,y2,nodes,edges2qmat,qvec)
x1 = -5.684;
y1 = -2.79;
x2 = -3.232;
y2 = -4.8;
epsilon = 0.002;
X = nodes(:,1);
Y = nodes(:,2);
a = find(X<x1+epsilon);
b = find(X>x1-epsilon);
xindices1 = intersect(a,b);
c = find(Y<y1+epsilon);
d = find(Y>y1-epsilon);
yindices1 = intersect(c,d);
index1 = intersect(xindices1,yindices1)
% f = find(X<x2+epsilon);
% g = find(X>x2-epsilon);
% xindices2 = intersect(f,g);
% h = find(Y<y2+epsilon);
% i = find(Y>y2-epsilon);
% yindices2 = intersect(h,i);
% index2 = intersect(xindices2,yindices2);
% e = edges2qmat(min(index1,index2),max(index1,index2));
% qvec(e)