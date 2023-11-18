function [Q,U,W,C,X,totiter] = pfqn_gflinearizer(L,N,Z,type,tol,maxiter,alpha)
R = size(L,2);
[Q,U,W,C,X,totiter] = pfqn_egflinearizer(L,N,Z,type,tol,maxiter,alpha*ones(1,R));
end


