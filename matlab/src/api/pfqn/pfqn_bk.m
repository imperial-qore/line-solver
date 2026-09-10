%{
%{
 % @file pfqn_bk.m
 % @brief Birman-Kogan saddle point normalizing constant with bottleneck detection.
%}
%}

%{
%{
 % @brief Birman-Kogan saddle point normalizing constant with bottleneck detection.
 %
 % Birman and Kogan (Stochastic Models 8(3):543-563, 1992) evaluate the
 % multichain partition function by the saddle point method applied to the
 % Cauchy inversion of its generating function. Stations that serve a single
 % chain and appear only once (the paper's dedicated single servers) stay
 % OUTSIDE the exponent as O(1) algebraic factors, so their poles may be
 % crossed by the saddle point; Algorithm 1 detects those chains and pins
 % their coordinate on the pole, which is where the residue rather than the
 % saddle carries the mass. The remaining stations are the paper's large
 % groups of identical stations and are exponentiated.
 %
 % @fn pfqn_bk(L, N, Z)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector (1 x classes).
 % @param Z Think time vector (default: zeros).
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
 % @return X Chain throughputs (the saddle point coordinates).
 % @return U Utilizations (stations x classes).
 % @return A Chains whose dedicated station is not saturated (eq. 29).
 % @return B Chains whose dedicated station is a bottleneck (eq. 30).
%}
%}
function [G,lG,X,U,A,B] = pfqn_bk(L,N,Z)
% [G,LG,X,U,A,B] = PFQN_BK(L,N,Z)
%
% Birman-Kogan saddle point expansion of the normalizing constant.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(Z)
    Z = zeros(1,size(L,2));
end
Z = sum(Z,1);
[M,R] = size(L);
if isempty(L) || isempty(N) || sum(N) == 0
    G = 1; lG = 0; X = zeros(1,R); U = zeros(M,R); A = 1:R; B = [];
    return
end
% An empty class contributes a factor of 1 and has no saddle coordinate: its
% own term -N_r*log(z_r) is absent, so the minimum runs to z_r = 0.
if any(N==0) && any(N>0)
    keep = N>0;
    [G,lG,Xk,Uk,Ak,Bk] = pfqn_bk(L(:,keep),N(keep),Z(keep));
    X = zeros(1,R); U = zeros(M,R);
    idx = find(keep);
    X(idx) = Xk; U(:,idx) = Uk; A = idx(Ak); B = idx(Bk);
    return
end
% Stations with no demand at all do not enter the generating function.
Lq = L;
Lq(sum(Lq,2) <= 0, :) = [];
Mq = size(Lq,1);

% Dedicated station of each chain: a station serving that chain alone and
% carrying no identical twin. A chain with several of them keeps the slowest,
% which is the pole the saddle point meets first.
%
% Three conditions guard the split, and each states one hypothesis of the
% paper rather than a numerical safeguard. The model must hold a group of
% identical stations, since it is only against M_j >> 1 replicas that a lone
% station is an O(1) factor instead of part of the exponent. The chain must
% have no think time, since the paper replaces the IS station with the
% dedicated servers and never has both: with a think time the exponent
% already confines the saddle and no pole can be crossed. And the station
% must serve that chain alone, which is what makes its pole lie on a
% coordinate axis, where Algorithm 1 can pin it.
mu = inf(1,R);
poleRow = zeros(1,R);
isPole = false(Mq,1);
mult = local_multiplicity(Lq);
if any(mult > 1)
    for ist = 1:Mq
        nz = find(Lq(ist,:) > 0);
        if numel(nz) ~= 1 || mult(ist) > 1
            continue
        end
        r = nz;
        if Z(r) > GlobalConstants.FineTol
            continue
        end
        if 1/Lq(ist,r) < mu(r)
            mu(r) = 1/Lq(ist,r);
            poleRow(r) = ist;
        end
    end
end
isPole(poleRow(poleRow>0)) = true;
Lg = Lq(~isPole,:); % the exponentiated stations (the paper's groups)

% Algorithm 1 as an active set method on the strictly convex psi. The bound
% z_r <= mu_r is the pole of chain r's dedicated station; a chain sitting on
% its bound is the paper's set B, for which the residue dominates.
[z,onBound] = local_minimize(Lg,N,Z,mu);
A = find(~onBound);
B = find(onBound);
X = z;
U = L .* repmat(X,M,1);
for r = B
    U(:,r) = min(U(:,r),1); % the dedicated station of a chain in B is saturated
end

psi0 = local_psi(z,Lg,N,Z);
if isempty(A) % eq. (25): every chain is pinned, the residues carry everything
    lG = psi0;
else
    H = local_hessian(z,Lg,N);
    HAA = H(A,A);
    [Hc,fail] = chol(HAA);
    if fail == 0
        logdetH = 2*sum(log(diag(Hc)));
    else
        [~,Uh,~] = lu(HAA);
        logdetH = sum(log(abs(diag(Uh))));
    end
    % eq. (22) and (24): the free coordinates carry the Gaussian prefactor and
    % the algebraic factor of their own dedicated station, the pinned ones do not
    lG = psi0 - 0.5*numel(A)*log(2*pi) - 0.5*logdetH - sum(log(z(A)));
    for r = A
        if isfinite(mu(r))
            lG = lG - log(1 - z(r)/mu(r));
        end
    end
end
G = exp(lG);
end

function mult = local_multiplicity(L)
% number of stations sharing each demand row, up to relative rounding
M = size(L,1);
mult = ones(M,1);
[Ls,ord] = sortrows(L);
s = 1;
for i = 2:M+1
    same = false;
    if i <= M
        scale = max([1, max(abs(Ls(i,:))), max(abs(Ls(i-1,:)))]);
        same = max(abs(Ls(i,:)-Ls(i-1,:))) <= GlobalConstants.FineTol*scale;
    end
    if ~same
        mult(ord(s:i-1)) = i-s;
        s = i;
    end
end
end

function f = local_psi(z,Lg,N,Z)
% exponent of the integrand, groups only (eq. 11 and 23 in unscaled variables)
f = Z*z' - N*log(z)';
if ~isempty(Lg)
    f = f - sum(log(1 - Lg*z'));
end
end

function g = local_grad(z,Lg,N,Z)
g = Z - N./z;
if ~isempty(Lg)
    g = g + (1./(1 - Lg*z'))'*Lg;
end
end

function H = local_hessian(z,Lg,N)
R = numel(z);
H = diag(N./z.^2);
if ~isempty(Lg)
    d = 1./(1 - Lg*z');
    H = H + Lg'*(repmat(d.^2,1,R).*Lg);
end
end

function [z,onBound] = local_minimize(Lg,N,Z,mu)
% minimize psi over {z>0, Lg*z'<1, z<=mu} by damped Newton on the free set
R = numel(N);
onBound = false(1,R);
z = local_init(Lg,N,Z,mu);
for outer = 1:(R+1)
    z(onBound) = mu(onBound);
    free = find(~onBound);
    if isempty(free)
        break
    end
    for it = 1:500
        g = local_grad(z,Lg,N,Z);
        if norm(g(free)) <= 1e-12*max(1,sum(N))
            break
        end
        H = local_hessian(z,Lg,N);
        dz = zeros(1,R);
        dz(free) = -(H(free,free)\g(free)')';
        alpha = 1;
        while true
            zt = z; zt(free) = z(free) + alpha*dz(free);
            if all(zt(free) > 0) && all(zt(free) <= mu(free)) && ...
                    (isempty(Lg) || max(Lg*zt') < 1)
                break
            end
            alpha = alpha/2;
            if alpha < 1e-14
                break
            end
        end
        if alpha < 1e-14
            break
        end
        z(free) = z(free) + alpha*dz(free);
    end
    % a chain whose descent direction still pushes past its pole belongs to B
    g = local_grad(z,Lg,N,Z);
    newly = free(z(free) >= mu(free)*(1-1e-9) & g(free) < 0);
    if isempty(newly)
        break
    end
    onBound(newly) = true;
end
z(onBound) = mu(onBound);
end

function z = local_init(Lg,N,Z,mu)
% start from the light traffic throughput, then retract into the domain
R = numel(N);
den = Z;
if ~isempty(Lg)
    den = den + sum(Lg,1);
end
z = N./max(den,GlobalConstants.FineTol);
z = min(z,0.99*mu);
if ~isempty(Lg)
    for it = 1:200
        if max(Lg*z') < 0.9
            break
        end
        z = z*0.7;
    end
end
end
