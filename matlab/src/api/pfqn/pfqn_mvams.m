%{
%{
 % @file pfqn_mvams.m
 % @brief General-purpose MVA for mixed networks with multiserver nodes.
%}
%}

%{
%{
 % @brief General-purpose MVA for mixed networks with multiserver nodes.
 % @fn pfqn_mvams(lambda, L, N, Z, mi, S)
 % @param lambda Arrival rate vector.
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mi Queue replication factors (default: ones).
 % @param S Number of servers per station (default: ones).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return CN Cycle times.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [XN,QN,UN,CN,lG]=pfqn_mvams(lambda,L,N,Z,mi,S)
% [XN,QN,UN,CN,LOGG]=PFQN_MVAMS(LAMBDA,L,N,Z,MI,S)

% this is a general purpose script to handle mixed qns with multi-server nodes
% S(i) number of servers in station i
[M,R]=size(L); % get number of queues (M) and classes (R)
Ntot = sum(N(isfinite(N)));
mu = ones(M,Ntot);
if nargin<6 %~exist('S','var')
    S = ones(M,1);
end
if nargin<5 %~exist('mi','var')
    mi = ones(M,1);
end
if isempty(Z)
    Z = zeros(1,R);
end
for ist=1:M
    mu(ist,:) = min(1:Ntot,S(ist)*ones(1,Ntot));
end
if max(S(isfinite(S))) == 1 % if no multi-server nodes
    if any(isinf(N)) % open or mixed model
        [XN,QN,UN,CN,lG] = pfqn_mvamx(lambda,L,N,Z,mi);
    else % closed model
        [XN,QN,UN,CN,lG] = pfqn_mva(L,N,Z,mi);
    end
else % if the model has multi-server nodes
    if any(isinf(N)) % open or mixed model
        if max(mi) == 1
            lG = NaN; % NC not available in this case
            [XN,QN,UN,CN] = pfqn_mvaldms(lambda,L,N,Z,S);
        else
            line_error(mfilename,'Queue replicas not available in exact MVA for mixed models.');
        end
    else
        [XN,QN,UN,CN,lG] = pfqn_mvald(L,N,Z,mu);
        lG=lG(end);
    end
end
return
end
