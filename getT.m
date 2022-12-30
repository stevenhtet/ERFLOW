function T = getT(factor)

lambda = 0.072527;
T = exprnd(lambda)/2;
while T<0.007 || T>0.5
    T = exprnd(lambda)/2;
end
T = T/factor;
% T = round(T,3)