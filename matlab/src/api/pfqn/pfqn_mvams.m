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
 % @return CN Residence times (M x R), as in PFQN_MVA.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [XN,QN,UN,CN,lG]=pfqn_mvams(lambda,L,N,Z,mi,S)
% [XN,QN,UN,CN,LOGG]=PFQN_MVAMS(LAMBDA,L,N,Z,MI,S)

% this is a general purpose script to handle mixed qns with multi-server nodes
% S(i) number of servers in station i
[M,R]=size(L); % get number of queues (M) and classes (R)
Ntot = 0;
hasOpenClasses = false;
for r = 1:R
    if isinf(N(r))
        hasOpenClasses = true;
    else
        Ntot = Ntot + N(r);
    end
end
mu = ones(M,Ntot);
if nargin<6 %~exist('S','var')
    S = ones(M,1);
end
if nargin<5 %~exist('mi','var')
    mi = ones(M,1);
end
if isempty(Z)
    Z = zeros(1,R);
elseif size(Z,1) > 1
    % Sum think times across multiple delay stations
    Z = sum(Z,1);
end
for ist=1:M
    mu(ist,:) = min(1:Ntot,S(ist)*ones(1,Ntot));
end
hasMultiServer = false;
for ist = 1:M
    if isfinite(S(ist)) && S(ist) > 1
        hasMultiServer = true;
        break;
    end
end

if ~hasMultiServer % if no multi-server nodes
    if hasOpenClasses % open or mixed model
        [XN,QN,UN,CN,lG] = pfqn_mvamx(lambda,L,N,Z,mi);
    else % closed model
        [XN,QN,UN,CN,lG] = pfqn_mva(L,N,Z,mi);
    end
else % if the model has multi-server nodes
    if hasOpenClasses % open or mixed model
        if max(mi) == 1
            lG = NaN; % NC not available in this case
            [XN,QN,UN,CN] = pfqn_mvaldms(lambda,L,N,Z,S);
        else
            line_error(mfilename,'Queue replicas not available in exact MVA for mixed models.');
        end
    else
        [XN,QN,UN,CN,lG] = pfqn_mvald(L,N,Z,mu);
        lG=lG(end);
        % PFQN_MVALD belongs to the load-dependent family and reports a (1 x R)
        % cycle time, but PFQN_MVAMS follows the PFQN_MVA contract as its other
        % three branches do: an (M x R) per-station residence time.
        CN = zeros(M,R);
        for r=1:R
            if N(r) > 0
                CN(:,r) = QN(:,r)/XN(r);
            else
                % an absent class, as in PFQN_MVA
                CN(:,r) = L(:,r).*mi(:);
            end
        end
    end
end
return
end
