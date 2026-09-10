%{
%{
 % @file pfqn_qsa.m
 % @brief Queue-Shift Approximation (QSA) for closed product-form networks.
%}
%}

%{
%{
 % @brief Queue-Shift Approximation (QSA) for closed product-form networks.
 %
 % Schweitzer, Serazzi and Broglia (Tools'98, LNCS 1469, pp. 267-279)
 % approximate the arrival-instant queue lengths through the absolute shift
 % Y_ri(K) = 1 + Q_i(K-e_r) - Q_i(K) of the aggregate queue length, in place of
 % the fractional deviations of Linearizer. The aggregate core problem (eq. 13a)
 % is imposed at K, at every K-e_s and, in the three-level variant of eq. (16),
 % at every K-e_s-e_t with the affine extrapolation of eq. (15).
 %
 % The quintuple (16) is solved as a single system by the damped Newton method
 % of Sect. 4. Successive substitution over the decomposed core problems is not
 % used: it is unstable near saturation and converges to the degenerate root in
 % which the bottleneck absorbs the whole population.
 %
 % @fn pfqn_qsa(L, N, Z, type, tol, maxiter, levels)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param type Scheduling strategy per station; SchedStrategy.INF marks a delay centre.
 % @param tol Residual tolerance (default: 1e-10).
 % @param maxiter Maximum number of Newton iterations (default: 100).
 % @param levels 2 for the two-level QSA of eq. (14), 3 for eq. (16) (default: 3).
 % @param QN0 (M x R) queue lengths that warm-start the Bard-Schweitzer initialization; empty for the default cold start.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return W Residence times.
 % @return C Cycle times.
 % @return X System throughput.
 % @return totiter Newton iterations performed.
%}
%}
function [Q,U,W,C,X,totiter] = pfqn_qsa(L,N,Z,type,tol,maxiter,levels,QN0)

[M,R] = size(L);
N = reshape(N,1,[]);
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
if nargin<4 || isempty(type)
    type = SchedStrategy.PS*ones(M,1);
end
if nargin<5 || isempty(tol)
    tol = 1e-10;
end
if nargin<6 || isempty(maxiter)
    maxiter = 100;
end
if nargin<7 || isempty(levels)
    levels = 3;
end
if nargin<8
    QN0 = [];
end
isQC = reshape(type,1,[]) ~= SchedStrategy.INF;
isQC = reshape(isQC(1:M),1,M);

Q = zeros(M,R);
U = zeros(M,R);
W = zeros(M,R);
C = zeros(1,R);
X = zeros(1,R);
totiter = 0;

if isempty(L) || all(max(L,[],1)==0) || all(N<=0)
    for r=1:R
        if N(r)>0 && Z(r)>0
            X(r) = N(r)/Z(r);
        end
        U(:,r) = X(r)*L(:,r);
    end
    return
end

% Populations touched by (16): K, every K-e_s, every K-e_s-e_t.
pops = N;
sIdx = zeros(1,R);
pIdx = zeros(R,R);
for s=1:R
    n = N; n(s) = n(s)-1;
    if all(n>=0)
        pops(end+1,:) = n; %#ok<AGROW>
        sIdx(s) = size(pops,1);
    end
end
if levels>=3
    for s=1:R
        for t=s:R
            n = N; n(s) = n(s)-1; n(t) = n(t)-1;
            if all(n>=0)
                pops(end+1,:) = n; %#ok<AGROW>
                pIdx(s,t) = size(pops,1);
                pIdx(t,s) = pIdx(s,t);
            end
        end
    end
end
nP = size(pops,1);

% Bard-Schweitzer at every population supplies the Newton starting point.
q = zeros(M,nP);
for p=1:nP
    q(:,p) = qsa_aggbs(L,pops(p,:),Z,isQC,QN0);
end
qc = find(isQC);
Ldc = zeros(1,R);
for r=1:R
    Ldc(r) = sum(L(~isQC,r));
end

x = reshape(q(qc,:),[],1);
n = numel(x);
[F,adm] = qsa_resid(x,L,Z,pops,sIdx,pIdx,qc,Ldc,levels);
fnrm = norm(F);
for totiter=1:maxiter
    if fnrm < tol
        break
    end
    J = zeros(n,n);
    for k=1:n
        h = 1e-7*max(1,abs(x(k)));
        xp = x; xp(k) = xp(k)+h;
        Fp = qsa_resid(xp,L,Z,pops,sIdx,pIdx,qc,Ldc,levels);
        J(:,k) = (Fp-F)/h;
    end
    lastwarn('');
    dx = -J\F;
    if any(~isfinite(dx))
        dx = -pinv(J)*F;
    end
    accepted = false;
    lambda = 1;
    for ls=1:40
        xn = x + lambda*dx;
        [Fn,admn] = qsa_resid(xn,L,Z,pops,sIdx,pIdx,qc,Ldc,levels);
        if admn && norm(Fn) < fnrm
            x = xn; F = Fn; fnrm = norm(Fn); adm = admn;
            accepted = true;
            break
        end
        lambda = lambda/2;
    end
    if ~accepted
        break
    end
end
if ~adm
    line_warning(mfilename,'QSA left the admissible region; the returned solution may be unphysical.');
end

% Disaggregate (13) at K into the per-class measures
q(qc,:) = reshape(x,numel(qc),nP);
Y0 = qsa_shift(q,pops,sIdx,pIdx,0,0,levels);
for r=1:R
    if N(r)<1
        continue
    end
    W(isQC,r) = L(isQC,r).*(q(isQC,1) + Y0(isQC,r));
    W(~isQC,r) = L(~isQC,r);
    X(r) = N(r)/(Z(r) + sum(W(:,r)));
    Q(:,r) = X(r)*W(:,r);
    U(:,r) = X(r)*L(:,r);
    C(r) = N(r)/X(r) - Z(r);
end
end

% ---------------------------------------------------------------------------
function [F,adm] = qsa_resid(x,L,Z,pops,sIdx,pIdx,qc,Ldc,levels)
% Residual of the core equations (13) imposed simultaneously at every
% population of (16). adm flags the side conditions of Remark 2 (non-negative
% queue lengths, positive cycle times).
[M,R] = size(L);
nP = size(pops,1);
mq = numel(qc);
q = zeros(M,nP);
q(qc,:) = reshape(x,mq,nP);
F = zeros(mq,nP);
adm = all(x>=0);
for p=1:nP
    np = pops(p,:);
    [s,t] = qsa_which(p,sIdx,pIdx);
    Y = qsa_shift(q,pops,sIdx,pIdx,s,t,levels);
    A = q(qc,p) + Y(qc,:);           % 1 + Q_i(K - e_r) at the arrival instant
    acc = zeros(mq,1);
    for r=1:R
        if np(r)<1
            continue
        end
        c = Z(r) + sum(L(qc,r).*A(:,r)) + Ldc(r);
        if ~(c>0) || ~isfinite(c)
            adm = false;
            c = eps;
        end
        acc = acc + (np(r)/c)*L(qc,r).*A(:,r);
    end
    F(:,p) = q(qc,p) - acc;
end
F = reshape(F,[],1);
end

% ---------------------------------------------------------------------------
function Y = qsa_shift(q,pops,sIdx,pIdx,s,t,levels)
% Shift matrix (M x R) of eq. (16d)-(16e) at population K (s=t=0), at K-e_s
% (t=0) or, via the extrapolation (15), at K-e_s-e_t.
M = size(q,1);
R = size(pops,2);
Y = zeros(M,R);
if s==0
    for r=1:R
        if sIdx(r)>0
            Y(:,r) = 1 + q(:,sIdx(r)) - q(:,1);
        end
    end
elseif t==0
    if levels<3
        Y = qsa_shift(q,pops,sIdx,pIdx,0,0,levels); % (14): Y(K-e_s) ~= Y(K)
        return
    end
    for r=1:R
        if pIdx(s,r)>0 && pops(sIdx(s),r)>=1
            Y(:,r) = 1 + q(:,pIdx(s,r)) - q(:,sIdx(s));
        end
    end
else
    Ys = qsa_shift(q,pops,sIdx,pIdx,s,0,levels);
    Yt = qsa_shift(q,pops,sIdx,pIdx,t,0,levels);
    Y0 = qsa_shift(q,pops,sIdx,pIdx,0,0,levels);
    Y = Ys + Yt - Y0;
end
end

% ---------------------------------------------------------------------------
function [s,t] = qsa_which(p,sIdx,pIdx)
% Decode population index p into the removed classes: (0,0) for K, (s,0) for
% K-e_s, (s,t) for K-e_s-e_t.
s = 0; t = 0;
if p==1
    return
end
k = find(sIdx==p,1);
if ~isempty(k)
    s = k;
    return
end
[ss,tt] = find(pIdx==p);
s = ss(1); t = tt(1);
end

% ---------------------------------------------------------------------------
function q = qsa_aggbs(L,n,Z,isQC,QN0)
% Aggregate queue lengths of the Bard-Schweitzer fixed point at population n,
% with the delay-centre demands folded into the think time.
[M,R] = size(L);
q = zeros(M,1);
n(n<0) = 0;
if all(n<=0)
    return
end
Zeff = Z;
for r=1:R
    Zeff(r) = Z(r) + sum(L(~isQC,r));
end
if any(isQC)
    if isempty(QN0)
        [X,QN] = pfqn_bs(L(isQC,:),n,Zeff);
    else
        [X,QN] = pfqn_bs(L(isQC,:),n,Zeff,1e-6,1000,QN0(isQC,:));
    end
    q(isQC) = sum(QN,2);
else
    X = zeros(1,R);
    for r=1:R
        if n(r)>=1 && Zeff(r)>0
            X(r) = n(r)/Zeff(r);
        end
    end
end
for i=reshape(find(~isQC),1,[])
    q(i) = sum(reshape(X,1,[]).*L(i,:));
end
end
