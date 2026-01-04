function [pi,pis,pi0,scc,isrec] = ctmc_solve_reducible(Q,pi0,options)
if nargin==1
    [pi,pis,pi0,scc,isrec] = dtmc_solve_reducible(ctmc_randomization(Q), [], struct("tol",1e-12));
elseif nargin==2
    [pi,pis,pi0,scc,isrec] = dtmc_solve_reducible(ctmc_randomization(Q), pi0, struct("tol",1e-12));
else
    [pi,pis,pi0,scc,isrec] = dtmc_solve_reducible(ctmc_randomization(Q), pi0, options);
end
end