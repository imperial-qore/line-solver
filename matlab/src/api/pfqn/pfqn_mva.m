%{
%{
 % @file pfqn_mva.m
 % @brief Exact Mean Value Analysis (MVA) for product-form queueing networks.
%}
%}

function [XN,QN,UN,CN,lGN] = pfqn_mva(L,N,Z,mi)
%{
%{
 % @brief Exact Mean Value Analysis (MVA) for product-form queueing networks.
 % @fn pfqn_mva(L, N, Z, mi)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R).
 % @param mi (Optional) Additive term of the residence-time recursion
 %        C(i,s)=L(i,s)*(mi(i)+Qarv), 1 for a queueing station. Default: ones.
 %        THIS IS NOT A SERVER COUNT: mi(i)=c INFLATES the residence time by c
 %        rather than adding c servers. For multiserver stations call
 %        PFQN_MVAMS(lambda,L,N,Z,mi,S), which passes S to the load-dependent
 %        recursion with mu(i,n)=min(n,S(i)).
 % @return XN System throughput (1 x R).
 % @return QN Mean queue length (M x R).
 % @return UN Utilization (M x R).
 % @return CN Residence time (M x R).
 % @return lGN Logarithm of the normalizing constant.
%}
%}
% [XN,QN,UN,CN,LGN] = PFQN_MVA(L,N,Z,MI)
% [XN,QN,UN,CN] = pfqn_mva(L,N,Z,mi)
%
% Standard arrival theorem. For the interlocked-flow correction of Franks
% (1999), Ch. 4, Eq. (4.7), call PFQN_MVA_ILOCK instead.
XN=[];
QN=[];
UN=[];
CN=[];
lGN = 0;
InfServ=1;
if nargin == 2
    InfServ=0;
end
N = ceil(N);
[M_original,R]=size(L); % M stations, R classes
N=N(:)';
if nargin<4
    mi=ones(1,M_original);
end
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end

% Station consolidation disabled: pfqn_unique merges stations with identical
% demand rows, but this is incorrect for tandem (serial) networks where distinct
% stations happen to have the same service demand. The consolidation treats them
% as replicated (parallel) copies, producing wrong queue lengths and response times.
M = M_original;
mapping = 1:M_original;
if (~any(N>0))
    %line_warning(mfilename,'closed populations are empty');
    return
end
NR=length(N);
if (R~=NR)
    line_error(mfilename,'demand matrix and population vector have different number of classes');
end

XN=zeros(1,R);
QN=zeros(M,R);
UN=zeros(M,R);
CN=zeros(M,R);
if InfServ==1
    Z=Z(:)';
else
    Z=zeros(1,R);
end

prods=zeros(1,R-1); % generate population indices
for w=1:R-1
    prods(1,w) = prod(ones(1,R-(w+1)+1)+N(1,w+1:R));
end

firstnonempty=R;
while (N(firstnonempty)==0)
    firstnonempty = firstnonempty-1;
end

totpop=prod(N+1);
ctr=totpop;
Q=zeros(totpop,M);
currentpop=2;

n=zeros(1,R);
n(1,firstnonempty)=1;
while ctr % for each population
    s=1;
    while s <= R
        pos_n_1s=0;
        if n(s)>0
            n(s) = n(s)-1;
            pos_n_1s= n(R);
            w=1;
            while w <= R-1
                pos_n_1s = pos_n_1s + n(w)*prods(w);
                w=w+1;
            end % while w <= R-1
            n(s) = n(s)+1;
        end % if
        CNtot=0;
        i=1;
        while i <= M
            Lis=L(i,s);
            qarv=Q(1+pos_n_1s,i);
            CN(i,s)=Lis*(mi(i)+qarv);
            CNtot=CNtot+CN(i,s);
            i=i+1;
        end % while i <= M
        XN(s)=n(s)/(Z(s)+CNtot);
        i=1;
        while i <= M
            QN(i,s)=XN(s)*CN(i,s);
            Q(currentpop,i)=Q(currentpop,i)+QN(i,s);
            i=i+1;
        end % while i <= M
        s=s+1;
    end % while s <= R
    s=R;
    while s>0 && (n(1,s)==N(s)) || s>firstnonempty
        s=s-1;
    end
    % now compute the normalizing constant
    last_nnz = last_nonzero_index(n);
    if last_nnz > 0 && sum(n(1:last_nnz-1)) == sum(N(1:last_nnz-1)) && sum(n((last_nnz+1):R))==0
        logX = log(XN(last_nnz));
        lGN = lGN - logX;
    end
    if s==0
        break;
    end
    n(s)=n(s)+1;
    s=s+1;
    while s<=R
        n(s)=0;
        s=s+1;
    end
    ctr=ctr-1;
    currentpop=currentpop+1;
end
for m=1:M
    for r=1:R
        UN(m,r)=XN(r)*L(m,r);
    end
end

% Expand results back to original dimensions if stations were consolidated
if M < M_original
    [QN, UN, CN] = pfqn_expand(QN, UN, CN, mapping, M_original);
end
end

function idx = last_nonzero_index(v)
idx = 0;
for i = length(v):-1:1
    if v(i) ~= 0
        idx = i;
        return
    end
end
end
