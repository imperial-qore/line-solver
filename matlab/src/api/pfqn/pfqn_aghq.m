%{
%{
 % @file pfqn_aghq.m
 % @brief Adaptive Gauss-Hermite quadrature for the normalizing constant.
%}
%}

%{
%{
 % @brief Adaptive Gauss-Hermite quadrature for the normalizing constant.
 % @fn pfqn_aghq(L, N, Z, q)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR).
 % @param q Nodes per simplex direction (default: 3). q=1 reproduces pfqn_le.
 % @return Gn Estimated normalizing constant.
 % @return lGn Logarithm of normalizing constant.
%}
%}
function [Gn,lGn]=pfqn_aghq(L,N,Z,q)
% [GN,LGN]=PFQN_AGHQ(L,N,Z,Q)

% PFQN_AGHQ Solution of closed product-form queueing networks by adaptive
% Gauss-Hermite quadrature of the McKenna-Mitra integral form.
%
% [Gn,lGn]=pfqn_aghq(L,N,Z,q)
% Input:
% L : MxR demand matrix. L(i,r) is the demand of class-r at queue i
% N : 1xR population vector. N(r) is the number of jobs in class r
% Z : 1xR think time vector. Z(r) is the total think time of class r
% q : number of nodes per simplex direction (default: 3)
%
% Output:
% Gn : estimated normalizing constant
% lGn: logarithm of Gn. If Gn exceeds the floating-point range, only lGn
%      will be correctly estimated.
%
% Method. At Z=0 the radius of the McKenna-Mitra integral separates exactly and
% leaves int_Delta prod_r (L_r'*x)^N_r dx over the simplex, which in logistic
% coordinates w is int exp(h(w)) dw with h(w)=sum_r N_r*log(L_r'*x(w))+sum_i
% log x_i. Rescaling by the mode w* and the curvature A of pfqn_le,
% w=w*+A^(-1/2)*z, and applying the q-node Gauss-Hermite rule of the
% probabilists' weight in each of the M-1 directions gives
%
%   int exp(h) dw = det(A)^(-1/2) * sum_k prod_i wt(k_i) *
%                   exp(h(w*+A^(-1/2)*z_k) + z_k'*z_k/2).
%
% At q=1 the single node is the mode and the weight sqrt(2*pi), so the rule
% collapses to pfqn_le, up to the tolerance of the shared fixed point: LE is the
% first term of a convergent quadrature rather than an approximation of unknown
% accuracy. The cost is q^(M-1) evaluations, which is what confines the method
% to small M.
%
% A tensor-product rule is not invariant to the choice of square root of A: any
% B with B*B'=inv(A) is admissible and they place the nodes differently. The
% principal-axis frame from eig is used here, as in the reference results; where
% two curvatures are close to equal the frame is close to arbitrary and two
% valid rules can part company well above their own error, converging back
% together as q grows. Compare across codebases at the level of the answer, not
% node by node.
%
% Convergence is fast when the integrand is close to Gaussian in w and slow in
% the corner regime, where the mode sits at O(1/N) from a vertex of the simplex
% and the non-bottleneck directions are locally Gamma(1) rather than Gaussian:
% there the constant of pfqn_ble is the better instrument.
%
% With Z>0 the integrand is no longer homogeneous, so the radius is integrated
% numerically instead and the rule is applied to the M-1 simplex directions of
% int_Delta J(L'*x) dx. Every node then costs one radial quadrature, so this
% branch is q^(M-1) times slower than the q=1 rule. Because the
% radius is integrated rather than Laplaced, q=1 there returns the logistic
% expansion with an exact radius, which is NOT pfqn_le's own Z>0 branch.
%
% References:
% J. McKenna, D. Mitra. Integral representations and asymptotic expansions for
% closed Markovian queueing networks: normal usage. Bell Syst. Tech. J. 61(5),
% 1982.
% G. Casale. Accelerating performance inference over closed systems by
% asymptotic methods. ACM SIGMETRICS 2017.

[M,R]=size(L);
N=N(:)';
if nargin<3 || isempty(Z)
    Z=zeros(1,numel(N));
else
    Z=sum(Z,1);
end
if nargin<4 || isempty(q)
    q=3;
end

if isempty(L) || isempty(N) || sum(N)==0 || sum(L(:))<1e-4
    lGn = -sum(factln(N)) + sum(N(N>0).*log(Z(N>0)));
    Gn = exp(lGn);
elseif sum(Z(:))<GlobalConstants.Zero
    umax = pfqn_aghq_fpi(L,N);
    A = pfqn_aghq_hessian(L,N,umax);
    ld = pfqn_aghq_logdet(A);
    S = 0;
    for r=1:R
        S = S + N(r)*log(umax'*L(:,r));
    end
    h0 = S + sum(log(umax));
    w0 = log(umax(1:M-1)/umax(M));
    lacc = pfqn_aghq_rule(@(w) pfqn_aghq_h(w,L,N),w0,h0,A,q,M-1);
    lGn = multinomialln([N,M-1]) + factln(M-1) + h0 + lacc - 0.5*ld;
    Gn = exp(lGn);
else % Z>0: radius integrated exactly, quadrature over the simplex
    [vg,wg] = pfqn_aghq_gausslegendre(64);
    [xmax,A,ld,h0] = pfqn_aghq_simplex(L,N,Z,vg,wg);
    w0 = log(xmax(1:M-1)/xmax(M));
    lacc = pfqn_aghq_rule(@(w) pfqn_aghq_hZ(w,L,N,Z,M,vg,wg),w0,h0,A,q,M-1);
    lGn = -sum(factln(N)) + h0 + lacc - 0.5*ld;
    Gn = exp(lGn);
end
end

function lacc=pfqn_aghq_rule(hfun,w0,h0,A,q,d)
% LACC=PFQN_AGHQ_RULE(HFUN,W0,H0,A,Q,D)

% Log of sum_k prod_i wt(k_i)*exp(h(w0+A^(-1/2)*z_k)-h0+z_k'*z_k/2) over the
% d-dimensional tensor grid, accumulated with a running maximum. The
% det(A)^(-1/2) of the rule is applied by the caller.
if d==0
    lacc = 0;
    return
end
nodes = q^d;
if nodes > 1e7
    line_error(mfilename,sprintf('the tensor rule needs q^(M-1)=%d nodes; reduce q or use pfqn_le.',nodes));
end
[V,D] = eig((A+A')/2);
lam = diag(D);
if min(lam) <= 0
    lacc = NaN;
    return
end
B = V*diag(1./sqrt(lam));
[z,wt] = pfqn_aghq_hermite(q);
lwt = log(wt);
idx = ones(d,1);
lmax = -Inf;
s = 0;
for k=1:nodes
    zz = z(idx);
    lt = sum(lwt(idx)) + hfun(w0+B*zz) - h0 + 0.5*(zz'*zz);
    if lt > lmax
        s = s*exp(lmax-lt) + 1;
        lmax = lt;
    else
        s = s + exp(lt-lmax);
    end
    for j=d:-1:1
        idx(j) = idx(j)+1;
        if idx(j) <= q
            break
        end
        idx(j) = 1;
    end
end
lacc = lmax + log(s);
end

function h=pfqn_aghq_h(w,L,N)
% H=PFQN_AGHQ_H(W,L,N)

% Log-integrand on the simplex in logistic coordinates, Z=0, Jacobian included.
x = pfqn_aghq_softmax(w);
h = sum(N.*log(x'*L)) + sum(log(x));
end

function h=pfqn_aghq_hZ(w,L,N,Z,M,vg,wg)
% H=PFQN_AGHQ_HZ(W,L,N,Z,M,VG,WG)

% Same with think times: the radial integral replaces the homogeneous factor.
x = pfqn_aghq_softmax(w);
h = pfqn_aghq_radial(x'*L,N,Z,M,vg,wg) + sum(log(x));
end

function x=pfqn_aghq_softmax(w)
% X=PFQN_AGHQ_SOFTMAX(W)

a = [w(:);0];
e = exp(a-max(a));
x = e/sum(e);
end

function [z,w]=pfqn_aghq_hermite(q)
% [Z,W]=PFQN_AGHQ_HERMITE(Q)

% Q-point Gauss-Hermite rule of the probabilists' weight exp(-z^2/2), by
% Golub-Welsch. The weights sum to sqrt(2*pi).
k = (1:q-1)';
[V,D] = eig(diag(sqrt(k),1)+diag(sqrt(k),-1));
[z,i] = sort(diag(D));
w = sqrt(2*pi)*(V(1,i).^2)';
end

function u=pfqn_aghq_fpi(L,N)
% U=PFQN_AGHQ_FPI(L,N)

% Mode of the simplex integrand in logistic coordinates, Z=0, as pfqn_le_fpi.
M = size(L,1);
eta = sum(N)+M;
u = ones(M,1)/M;
u_1 = Inf*u;
it = 0;
while norm(u-u_1,1)>1e-11 && it<100000
    u_1 = u;
    c = u_1'*L;
    u = (1+u_1.*(L*(N./c)'))/eta;
    it = it+1;
end
end

function A=pfqn_aghq_hessian(L,N,u)
% A=PFQN_AGHQ_HESSIAN(L,N,U)

% Reduced (M-1)x(M-1) curvature at the mode, as pfqn_le_hessian.
M = size(L,1);
eta = sum(N)+M;
c = u'*L;
B = L.*(sqrt(N)./c);
K = B*B';
s = L*(N./c)';
A = (u*u').*(K-eta);
A(1:M+1:end) = eta*u.*(1-u) - u.*s + (u.^2).*diag(K);
A = A(1:M-1,1:M-1);
end

function u=pfqn_aghq_fpiZ(L,N,Z)
% U=PFQN_AGHQ_FPIZ(L,N,Z)

% Logistic-expansion mode with think times, warm start for the exact-radius
% fixed point. As pfqn_le_fpiZ, with v=eta-sum_r xi_r*Z_r.
M = size(L,1);
eta = sum(N)+M;
u = ones(M,1)/M;
v = eta;
u_1 = Inf*u;
v_1 = Inf;
it = 0;
while norm(u-u_1,1)+abs(v-v_1)>1e-11 && it<100000
    u_1 = u;
    v_1 = v;
    c = Z + v*(u_1'*L);
    u = (1+u_1.*((Z+v*L)*(N./c)'))/eta;
    v = eta - sum((N./c).*Z);
    it = it+1;
end
end

function [x,A,ld,h0]=pfqn_aghq_simplex(L,N,Z,vg,wg)
% [X,A,LD,H0]=PFQN_AGHQ_SIMPLEX(L,N,Z,VG,WG)

% Mode and curvature of h(w)=log J(L'*x(w))+sum_i log x_i with J the exact
% radial integral; pfqn_aghq_fpiZ is the fixed point and pfqn_aghq_hessian the
% curvature identity.
M = size(L,1);
x = pfqn_aghq_fpiZ(L,N,Z);
x_1 = Inf*x;
it = 0;
while norm(x-x_1,1)>1e-11 && it<10000
    x_1 = x;
    [~,G,vbar] = pfqn_aghq_radial(x_1'*L,N,Z,M,vg,wg);
    x = (1+x_1.*(L*G'))/vbar;
    x = x/sum(x);
    it = it+1;
end
[lJ,~,~,Lam] = pfqn_aghq_radial(x'*L,N,Z,M,vg,wg);
P = L*Lam*L' - diag(1./x.^2);
Jm = diag(x) - x*x';
Jm = Jm(:,1:M-1);
A = -Jm'*P*Jm;
A = (A+A')/2;
h0 = lJ + sum(log(x));
ld = pfqn_aghq_logdet(A);
end

function [lJ,G,vbar,Lam]=pfqn_aghq_radial(c,N,Z,M,vg,wg)
% [LJ,G,VBAR,LAM]=PFQN_AGHQ_RADIAL(C,N,Z,M,VG,WG)

% log J(c)=log int_0^inf exp(-v)*v^(M-1)*prod_r (Z_r+v*c_r)^N_r dv and the
% moments of the tilted law of v.
R = numel(N);
t = log(sum(N)+M);
for it=1:200
    v = exp(t);
    d = max(Z+v*c, realmin);
    F1 = -v + M + sum(N.*(v*c)./d);
    F2 = -v + sum(N.*(v*c).*Z./d.^2);
    if F2 > -1e-300
        break
    end
    step = -F1/F2;
    step = max(min(step,2),-2);
    if abs(step) < 1e-13
        t = t+step;
        break
    end
    t = t+step;
end
v = exp(t);
d = max(Z+v*c, realmin);
F2 = -v + sum(N.*(v*c).*Z./d.^2);
if F2 < -1e-300
    sig = 1/sqrt(-F2);
else
    sig = 1;
end
Fm = pfqn_aghq_logf(t,c,N,Z,M);
a = min(12*sig, t+745);
for k=1:60
    if t-a <= -745 || pfqn_aghq_logf(t-a,c,N,Z,M) < Fm-60
        break
    end
    a = min(1.6*a, t+745);
end
b = 12*sig;
for k=1:60
    if pfqn_aghq_logf(t+b,c,N,Z,M) < Fm-60
        break
    end
    b = 1.6*b;
end
tt = [0.5*a*vg + (t-0.5*a); 0.5*b*vg + (t+0.5*b)];
W = [0.5*a*wg; 0.5*b*wg];
vv = exp(tt);
D = max(Z+vv*c, realmin);
Fv = -vv + M*tt + log(D)*N';
mx = max(Fv);
e = W.*exp(Fv-mx);
se = sum(e);
lJ = mx + log(se);
p = e/se;
T = (vv.*N)./D;
G = p'*T;
vbar = p'*vv;
ET2 = T'*(p.*T);
dg = zeros(R,1);
Nc = N(:);
nz = Nc>0;
d2 = diag(ET2);
dg(nz) = d2(nz)./Nc(nz);
Lam = ET2 - G'*G - diag(dg);
Lam = (Lam+Lam')/2;
end

function F=pfqn_aghq_logf(t,c,N,Z,M)
% F=PFQN_AGHQ_LOGF(T,C,N,Z,M)

% Log-integrand of the radial integral in t=log v, Jacobian included.
v = exp(t);
F = -v + M*t + sum(N.*log(max(Z+v*c,realmin)));
end

function [x,w]=pfqn_aghq_gausslegendre(n)
% [X,W]=PFQN_AGHQ_GAUSSLEGENDRE(N)

% N-point Gauss-Legendre rule on [-1,1] by Golub-Welsch.
k = (1:n-1)';
b = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(b,1)+diag(b,-1));
[x,i] = sort(diag(D));
w = 2*(V(1,i).^2)';
end

function ld=pfqn_aghq_logdet(A)
% LD=PFQN_AGHQ_LOGDET(A)

[C,p] = chol(A);
if p==0
    ld = 2*sum(log(diag(C)));
else
    ld = log(det(A));
end
end

function mln=multinomialln(n)
% MLN=MULTINOMIALLN(N)

mln = factln(sum(n))- sum(factln(n));
end

function lf=factln(n)
% LF=FACTLN(N)

lf = gammaln(1+n);
end
