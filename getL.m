function L = getL()
%draw pinch length from normal distribution with Holcman 2018 parameters
L = normrnd(0.12,0.0395)/2;
while L<0.025 || L>0.24 
L = normrnd(0.12,0.0395)/2;
end