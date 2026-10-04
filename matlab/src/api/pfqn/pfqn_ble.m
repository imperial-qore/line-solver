%{
%{
 % @file pfqn_ble.m
 % @brief Logistic expansion with the eps->0 bias correction (LE+).
%}
%}

%{
%{
 % @brief Logistic expansion with the eps->0 bias correction (LE+).
 % @fn pfqn_ble(L, N, Z)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR).
 % @return Gn Estimated normalizing constant.
 % @return lGn Logarithm of normalizing constant.
%}
%}
function [Gn,lGn]=pfqn_ble(L,N,Z)
% [GN,LGN]=PFQN_BLE(L,N,Z)
%
% PFQN_BLE Asymptotic solution of closed product-form queueing networks by
% logistic expansion, corrected by the deficit that the expansion carries on the
% one model for which its reference measure is exact.
%
% Cas17 Theorem 2 holds for eps >= eps_N > 0; the K(1+eps*N) self-looping
% populations are what make the integrand concentrate. Evaluated at eps->0, as
% pfqn_le does, the curvature at the saddle tends to 1 rather than growing with
% N, so Laplace's method has no asymptotic regime there and carries an O(1)
% relative bias per Laplaced direction. The size of that bias follows from the
% control-variate reading of the closure: closing an integral with a reference
% measure whose exact value is known corrects the Laplace value by the factor
% the reference itself misses, and for a Dirichlet reference that factor is
%
%   -log L_alpha = sum_k r(alpha_k) - r(alpha_0),
%   r(a) = gammaln(a) - ((a-1/2)*log(a) - a + log(2*pi)/2),
%
% the Stirling remainder of its own parameters, with kappa = r(1) = 1-log(2*pi)/2.
% The reference here is not fitted but read off the balanced model, where the
% integrand is constant on the simplex and the uniform measure alpha = 1 is
% exactly right. That gives the two branch constants:
%
%   Z=0 : M*kappa - r(M),  the M-1 simplex directions plus the Jacobian, with
%         the radial integral done exactly as Gamma(N+M);
%   Z>0 : M*kappa,         the same simplex term plus the radial direction,
%         which is Laplaced in t=log(v) and contributes exactly +r(M), so the
%         -r(M) cancels.
%
% Both are exact deficits on the balanced model rather than fitted constants.
% Measured against exact convolution over random Z=0 models the residual after
% M*kappa-r(M) has median +0.006 nats, against +0.067 for the (M-1)*kappa used
% previously on this branch. The published expansion is NOT in error; see
% _kb/03-api-layer.md.
%
% Input:
% L : MxR demand matrix. L(i,r) is the demand of class-r at queue i
% N : 1xR population vector. N(r) is the number of jobs in class r
% Z : 1xR think time vector. Z(r) is the total think time of class r
%
% Output:
% Gn : estimated normalizing constant
% lGn: logarithm of Gn. If Gn exceeds the floating-point range, only lGn
%      will be correctly estimated.
%
% Reference:
% G. Casale. Accelerating performance inference over closed systems by
% asymptotic methods. ACM SIGMETRICS 2017.

[M,R]=size(L);

if isempty(L) || isempty(N) || sum(N)==0 || sum(L(:))<1e-4
    % Z may be absent on this branch, and an empty class contributes 0, not 0*log(0).
    if nargin<3 || isempty(Z)
        Zt = zeros(1,numel(N));
    else
        Zt = sum(Z,1);
    end
    lGn = - sum(factln(N));
    for r=1:numel(N)
        if N(r)>0
            lGn = lGn + N(r)*log(Zt(r));
        end
    end
    Gn=exp(lGn);
elseif nargin<3 || isempty(Z) || sum(Z(:))<GlobalConstants.Zero
    umax=pfqn_ble_fpi(L,N);
    A=pfqn_ble_hessian(L,N,umax');
    S=0;
    for r=1:R
        S=S+N(r)*log(umax'*L(:,r));
    end
    % eps->0 bias correction: M*kappa-r(M), the exact deficit of pfqn_le on the
    % balanced model, where the uniform reference is correct. See the header and
    % _kb/03-api-layer.md. The radial integral is exact here as Gamma(N+M).
    lGn = multinomialln([N,M-1]) + factln(M-1) + (M-1)*log(sqrt(2*pi)) + M*(1-log(2*pi)/2) - stirlingrem(M) - log(sqrt(det(A))) + sum(log(umax)) + S;
    Gn=exp(lGn);
else % Z>0
    [umax,vmax]=pfqn_ble_fpiZ(L,N,Z);
    A=pfqn_ble_hessianZ(L,N,Z,umax',vmax);
    S=0;
    for r=1:R
        S=S+N(r)*log(Z(r)+vmax*umax'*L(:,r));
    end
    % M Laplaced directions here, the M-1 of the simplex plus the radius. The
    % radial one is Laplaced in t=log(v) and contributes exactly +r(M), which
    % cancels the -r(M) of the Z=0 branch above, leaving M*kappa. See the header.
    lGn = -sum(factln(N)) -vmax + M*log(vmax) + M*log(sqrt(2*pi)) + M*(1-log(2*pi)/2) - log(sqrt(det(A))) + sum(log(umax)) + S;
    Gn=exp(lGn);
end
end

function rem=stirlingrem(a)
% REM=STIRLINGREM(A)

% Stirling remainder r(a)=log(Gamma(a))-((a-1/2)*log(a)-a+log(2*pi)/2), the amount
% by which Laplace's method underestimates log(Gamma(a)). r(1)=1-log(2*pi)/2 and
% r(a)=1/(12*a)+O(a^-2).
rem = gammaln(a) - ((a-0.5).*log(a) - a + 0.5*log(2*pi));
end

function [u,d]=pfqn_ble_fpi(L,N)
% find location of mode of gaussian
[M,R]=size(L);
u=ones(M,1)/M;
u_1=Inf*u;
d=[];
while norm(u-u_1,1)>1e-10
    u_1=u;
    for i=1:M
        u(i)=1/(sum(N)+M);
        for r=1:R
            u(i)=u(i)+N(r)/(sum(N)+M)*L(i,r)*u_1(i)/(u_1'*L(:,r));
        end
    end
    d=[d; abs(u-u_1)']; %#ok<AGROW>
end
end

function [u,v,d]=pfqn_ble_fpiZ(L,N,Z)
% find location of mode of gaussian
[M,R]=size(L);
eta = sum(N)+M;
u=ones(M,1)/M;
% Note: eq. (35) in the SIGMETRICS 2017 paper has a spurious +1 in the v
% equation; the correct stationary point is v = eta - sum_r xi_r*Z_r.
v=eta;
u_1=Inf*u;
v_1=Inf*v; %#ok<NASGU>
d=[];
while norm(u-u_1,1)>1e-10
    u_1=u;
    v_1=v;
    for ist=1:M
        u(ist)=1/eta;
        for r=1:R
            u(ist)=u(ist)+(N(r)/eta)*(Z(r)+v*L(ist,r))*u_1(ist)/(Z(r)+v*u_1'*L(:,r));
        end
    end
    xi = zeros(1, R);
    for r=1:R
        xi(r)=N(r)/(Z(r)+v*u_1(:)'*L(:,r));
    end
    v=eta;
    for r=1:R
        v=v-xi(r)*Z(r);
    end
    d=[d; abs(u-u_1)'+abs(v-v_1)]; %#ok<AGROW>
end
end

function hu=pfqn_ble_hessian(L,N,u0)
% find hessian of gaussian
[M,R]=size(L);
Ntot=sum(N);
hu=zeros(M-1);
for i=1:(M-1)
    for j=1:(M-1)
        if i~=j
            hu(i,j)=-(Ntot+M)*u0(i)*u0(j);
            for r=1:R
                hu(i,j)=hu(i,j)+N(r)*L(i,r)*L(j,r)*(u0(i)*u0(j))/(u0*L(:,r))^2;
            end
        else
            hu(i,j)=(Ntot+M)*u0(i)*sum(allbut(u0,i));
            for r=1:R
                hu(i,j)=hu(i,j)-N(r)*L(i,r)*u0(i)*(allbut(u0,i)*L(allbut(1:M,i),r))/(u0*L(:,r))^2;
            end
        end
    end
end
end

function A=pfqn_ble_hessianZ(L,N,Z,u,v)
% find hessian of gaussian
[K,R]=size(L);
Ntot=sum(N);
A=zeros(K);
csi = zeros(1,R);
csi2N = zeros(1,R);
for r=1:R
    csi(r)=N(r)/(Z(r)+v*u*L(:,r));
    % csi(r)^2/N(r) rewritten as N(r)/c(r)^2. Identical where both are defined, but 0
    % rather than 0/0 for an empty class, which oner() makes routine in pfqn_nc.
    csi2N(r)=N(r)/(Z(r)+v*u*L(:,r))^2;
end
Lhat = zeros(K,R);
for k=1:K
    for r=1:R
        Lhat(k,r)=Z(r)+v*L(k,r);
    end
end
eta=Ntot+K;
for i=1:K
    for j=1:K
        if i~=j
            A(i,j)=-eta*u(i)*u(j);
            for r=1:R
                A(i,j)=A(i,j)+csi2N(r)*Lhat(i,r)*Lhat(j,r)*(u(i)*u(j));
            end
        end
    end
end
for i=1:K
    A(i,i)=-sum(allbut(A(i,:),i));
end
Ared=A(1:(K-1),1:(K-1));
A=zeros(K,K);
A(1:(K-1),1:(K-1))=Ared;
A(K,K)=1;
for r=1:R
    A(K,K)=A(K,K)-csi2N(r)*Z(r)*u*L(:,r);
end
A(K,K)=v*A(K,K);
for i=1:(K-1)
    A(i,K)=0;
    for r=1:R
        A(i,K)=A(i,K)+v*u(i)*(csi2N(r)*Lhat(i,r)*(u*L(:,r))-csi(r)*L(i,r));
    end
    A(K,i)=A(i,K);
end
end

function y=allbut(y,xset)
y=y(setdiff(1:length(y),xset));
end

function mln=multinomialln(n)
mln = factln(sum(n))- sum(factln(n));
end

function lf=factln(n)
lf = gammaln(1+n);
end
