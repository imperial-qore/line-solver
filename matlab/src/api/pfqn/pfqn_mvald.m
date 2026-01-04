%{
%{
 % @file pfqn_mvald.m
 % @brief Exact MVA for load-dependent closed queueing networks.
%}
%}

%{
%{
 % @brief Exact MVA for load-dependent closed queueing networks.
 % @fn pfqn_mvald(L, N, Z, mu, stabilize)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mu Load-dependent rate matrix (MxNt).
 % @param stabilize Force non-negative probabilities (default: true).
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return CN Cycle times.
 % @return lGN Logarithm of normalizing constant evolution.
 % @return isNumStable Numerical stability flag.
 % @return pi Marginal queue-length probabilities.
%}
%}
function [XN,QN,UN,CN,lGN,isNumStable,pi]=pfqn_mvald(L,N,Z,mu,stabilize)
% [XN,QN,UN,CN,LGN]=PFQN_MVALD(L,N,Z,MU)

[M,R]=size(L); % get number of queues (M) and classes (R)

XN = zeros(M,R);
QN = zeros(M,R);
UN = zeros(M,R);
CN = zeros(M,R);
WN = zeros(M,R);
lGN = -Inf;

if any(N<0)
    isNumStable = true;
    return
end
% stabilize ensures that probabilities do not become negative
if nargin<5
    stabilize = true;
end

warn = true;
isNumStable = true;
Xs=zeros(R,prod(N+1)); % throughput for a model with station i less
pi=ones(M,sum(N)+1,prod(N+1)); % marginal queue-length probabilities pi(k)
n=pprod(N); % initialize the current population
lGN = 0;
while n~=-1
    WN=0*WN;
    for s=1:R
        if n(s)>0
            for ist=1:M
                WN(ist,s)=0;
                for k=1:sum(n)
                    WN(ist,s)=WN(ist,s)+(L(ist,s)/mu(ist,k))*k*pi(ist,(k-1)+1,hashpop(oner(n,s),N));
                end
            end
            Xs(s,hashpop(n,N))=n(s)/(Z(s)+sum(WN(:,s)));
        end
    end
    
    % compute pi(k|n)
    for k=1:sum(n)
        for ist=1:M
            pi(ist,(k)+1,hashpop(n,N))=0;
        end
        for s=1:R
            if n(s)>0
                for ist=1:M
                    pi(ist,(k)+1,hashpop(n,N)) = pi(ist,(k)+1,hashpop(n,N)) + (L(ist,s)/mu(ist,k))*Xs(s,hashpop(n,N))*pi(ist,(k-1)+1,hashpop(oner(n,s),N));
                end
            end
        end
    end
        
    % compute pi(0|n)
    for ist=1:M
        p0 = 1-sum(pi(ist,(1:sum(n))+1,hashpop(n,N)));
        if p0<0
            if warn
                line_warning(mfilename,'MVA-LD is numerically unstable on this model, LINE will force all probabilities to be non-negative.\n'); 
%                N
                warn=false;
                isNumStable = false;
            end
            if stabilize
                pi(ist,(0)+1,hashpop(n,N)) = eps;
            else
                pi(ist,(0)+1,hashpop(n,N)) = p0;                
            end
        else
            pi(ist,(0)+1,hashpop(n,N)) = p0;
        end
    end
    
    last_nnz = find(n>0, 1, 'last' );
    if sum(n(1:last_nnz-1)) == sum(N(1:last_nnz-1)) && sum(n((last_nnz+1):R))==0
        logX = log(Xs(last_nnz,hashpop(n,N)));
        %hashpop(n,N)
        if ~isempty(logX)
            lGN(end+1) = lGN(end) - logX;
        end
    end
    
    n=pprod(n,N); % get the next population
end
X=Xs(:,hashpop(N,N));
XN=X';
pi=pi(:,:,hashpop(N,N));
QN = WN.*repmat(XN,M,1);
%UN = repmat(XN,M,1) .* L;
UN = 1-pi(:,1);
CN = N./XN - Z; % cycle time exclusive of think time
end

function i=hashpop(n,N)
% I=HASHPOP(N,N)

i=1; % index of the empty population
for r=1:length(N)
    i=i+prod(N(1:r-1)+1)*n(r);
end
end
