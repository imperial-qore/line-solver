%{
%{
 % @file pfqn_chow.m
 % @brief JMT-compatible Chow approximate MVA.
%}
%}

function [XN,QN,UN,RN,it]=pfqn_chow(L,N,Z,tol,maxiter,QN0,type,variant)
%{
%{
 % @brief JMT-compatible Chow approximate MVA.
 %
 % JMT's Chow analyzer estimates every class's arrival-instant queue length
 % at a station by the aggregate queue length at the full population. This is
 % the Bard large-customer-population fixed point implemented by pfqn_lcp.
 %
 % @fn pfqn_chow(L, N, Z, tol, maxiter, QN0, type, variant)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param tol Tolerance for convergence.
 % @param maxiter Maximum number of iterations.
 % @param QN0 Initial guess for queue lengths.
 % @param type Scheduling strategy type (default: PS).
 % @param variant Legacy compatibility argument; accepted but ignored.
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return it Number of iterations performed.
%}
%}

if nargin<3 || isempty(Z)
    Z=0*N;
end
if nargin<4 || isempty(tol)
    tol = 1e-6;
end
if nargin<5 || isempty(maxiter)
    maxiter = 1000;
end
[M,~]=size(L);
if nargin<6
    QN0 = [];
end
if nargin<7 || isempty(type)
    type = SchedStrategy.PS * ones(M,1);
end
if nargin<8
    variant = []; %#ok<NASGU>
end

[XN,QN,UN,RN,it] = pfqn_lcp(L,N,Z,tol,maxiter,QN0,type);
end
