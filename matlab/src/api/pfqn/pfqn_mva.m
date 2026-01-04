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
 % @param mi (Optional) Server multiplicity vector (1 x M). Default: single servers.
 % @return XN System throughput (1 x R).
 % @return QN Mean queue length (M x R).
 % @return UN Utilization (M x R).
 % @return CN Residence time (M x R).
 % @return lGN Logarithm of the normalizing constant.
%}
%}
% [XN,QN,UN,CN,LGN] = PFQN_MVA(L,N,Z,MI)
% [XN,QN,UN,CN] = pfqn_mva(L,N,Z,mi)
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

% Detect and consolidate replicated stations
[L, ~, ~, mi_unique, mapping] = pfqn_unique(L);
[M,~] = size(L);

% Combine user-provided mi with detected multiplicity
% For each unique station j, sum the mi values of all original stations mapping to it
mi_combined = zeros(1, M);
for i = 1:M_original
    mi_combined(mapping(i)) = mi_combined(mapping(i)) + mi(i);
end
mi = mi_combined;
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
            CN(i,s)=Lis*(mi(i)+Q(1+pos_n_1s,i));
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
    last_nnz = find(n>0, 1, 'last' );
    if sum(n(1:last_nnz-1)) == sum(N(1:last_nnz-1)) && sum(n((last_nnz+1):R))==0
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
