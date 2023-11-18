function value = softmin(x,y,alpha) 
% VALUE = SOFTMIN(X,Y,ALPHA)
value = - ((-x)*exp(-alpha*x) -y*exp(-alpha*y)) / (exp(-alpha*x) + exp(-alpha*y));
end