%{
%{
 % @file pfqn_mva_ilock.m
 % @brief Exact MVA recursion carrying the interlocked-flow correction.
%}
%}

function [XN,QN,UN,CN,lGN] = pfqn_mva_ilock(L,N,Z,mi,IL)
%{
%{
 % @brief Exact MVA recursion carrying the interlocked-flow correction.
 % @fn pfqn_mva_ilock(L, N, Z, mi, IL)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R).
 % @param mi (Optional) Server multiplicity vector (1 x M). Default: single servers.
 % @param IL Interlock matrix (R x R). IL(r,s) is the share of the class-s queue
 %        that a class-r arrival cannot see, because that work was itself caused
 %        by the class-r request (Franks 1999, Ch. 4, Eq. 4.7).
 % @return XN System throughput (1 x R).
 % @return QN Mean queue length (M x R).
 % @return UN Utilization (M x R).
 % @return CN Residence time (M x R).
 % @return lGN Always NaN: the interlock leaves the model outside product form.
%}
%}
% [XN,QN,UN,CN,LGN] = PFQN_MVA_ILOCK(L,N,Z,MI,IL)
%
% Closed single-server models only. The correction replaces the arrival theorem
% term Q(n-1_s,i) by a per-class weighted sum sum_r ILw(s,r)*Qc(n-1_s,i,r), so
% the recursion has to carry per-class queue lengths that PFQN_MVA does not
% need. Pass an empty IL to PFQN_MVA instead.
%
% The discounted arrival-instant queue is floored at the in-service component,
% as in lqns MVA::queueOnly_adjusted, so the correction damps itself out as a
% station saturates. That is a self-limiting guard, NOT a hard capacity test:
% sum_s XN(s)*L(i,s) <= mi(i) is still asserted nowhere.
% See git show 449847e7b:_kb/log.md.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

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
if nargin<4 || isempty(mi)
    mi=ones(1,M_original);
end
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
if nargin<5 || isempty(IL)
    line_error(mfilename,'an interlock matrix is required; use pfqn_mva for the standard arrival theorem');
end
if size(IL,1)~=R || size(IL,2)~=R
    line_error(mfilename,'the interlock matrix must be nclasses x nclasses');
end
ILw = max(0,min(1,1-IL)); % share of each class queue that an arrival still sees
ILw(1:(R+1):end) = 1; % a request always sees its own class in full

M = M_original;
mapping = 1:M_original;
if (~any(N>0))
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
Qc=zeros(totpop,M,R); % per-class queue lengths, needed by the interlock
Uc=zeros(totpop,M,R); % per-class in-service component, the interlock's floor
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
            qarv=0;
            r=1;
            while r <= R
                % In-service protection, as in lqns MVA::queueOnly_adjusted:
                % the discount bites on the WAITING part only, never on the job
                % already in service. As a station saturates, Uc approaches Qc
                % and the correction damps itself out, which is what keeps the
                % interlocked rates near the station's capacity.
                qarv=qarv+max(ILw(s,r)*Qc(1+pos_n_1s,i,r), Uc(1+pos_n_1s,i,r));
                r=r+1;
            end
            CN(i,s)=Lis*(mi(i)+qarv);
            CNtot=CNtot+CN(i,s);
            i=i+1;
        end % while i <= M
        XN(s)=n(s)/(Z(s)+CNtot);
        i=1;
        while i <= M
            QN(i,s)=XN(s)*CN(i,s);
            Q(currentpop,i)=Q(currentpop,i)+QN(i,s);
            Qc(currentpop,i,s)=QN(i,s);
            Uc(currentpop,i,s)=XN(s)*L(i,s);
            i=i+1;
        end % while i <= M
        s=s+1;
    end % while s <= R
    s=R;
    while s>0 && (n(1,s)==N(s)) || s>firstnonempty
        s=s-1;
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
lGN = NaN; % the interlock leaves the model outside product form

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
