%{
%{
 % @file pfqn_kt.m
 % @brief Knessl-Tier asymptotic expansion for normalizing constant.
%}
%}

%{
%{
 % @brief Knessl-Tier asymptotic expansion for normalizing constant.
 % @fn pfqn_kt(L, N, Z)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector (default: zeros).
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
 % @return X System throughput.
 % @return Q Mean queue lengths.
%}
%}
function [G,lG,X,Q]=pfqn_kt(L,N,Z)
% Knessl-Tier asymptotic expansion, fixed to include the IS/think-time term
% (defect G18) and to evaluate the exponent at the exact saddle point.
%
% Derivation. In LINE's convention the generating function of G over the
% population vector is
%   sum_N G(N) prod_r u_r^N_r = exp(sum_r Z_r u_r) prod_k (1-sum_r L_kr u_r)^-1
% so Cauchy extraction and steepest descent on
%   F(u) = sum_r Z_r u_r - sum_k log(1-U_k) - sum_r N_r log u_r,  U_k = L(k,:)*u
% give
%   G ~ (2 pi)^{-R/2} det(H)^{-1/2} exp(F(u*)) / prod_r u*_r
%   H_rs = delta_rs N_r/u_r^2 + sum_k L_kr L_ks/(1-U_k)^2
% where u* solves N_r = u_r Z_r + sum_k u_r L_kr/(1-U_k) (asymptotic MVA
% fixed point). This is Knessl-Tier Result 2 rewritten in LINE's convention:
% their rho_k absorbs the think rate and the factor exp(sum_r Z_r u_r) is
% contained in Psi(1,y*) via the M!/(M-n)! Stirling terms. Stock pfqn_kt
% dropped the linear think term (it is absent from both the exponent and the
% Hessian) and evaluated the exponent at the AQL throughput.
if isempty(L) || isempty(N) || sum(N)==0
    G = 1;
    lG = 0;
    X = [];
    Q = [];
    return;
end
if nargin<3
    Z = N*0;
end
[Morig,Rorig] = size(L); %#ok<ASGLU>
% fix self-looping customers as they would yield Uk=1
slcdemandfactor = 0;
if Rorig>1
    isslc = false(1,Rorig);
    for r=1:Rorig
        if nnz(L(:,r))==1
            if Z(r)==0
                ist = find(L(:,r)>0);
                L = [L; repmat(L(ist,:),N(r),1)];
                isslc(r) = true;
                slcdemandfactor = N(r)*log(L(ist,r));
            end
        end
    end
    L(:,isslc)=[];
    Z(:,isslc)=[];
    N(:,isslc)=[];
end
[M,R] = size(L);
Ntot = sum(N);
if Ntot <= 4
    [X,Q] = pfqn_bs(L,N,Z);
else
    [X,Q] = pfqn_aql(L,N,Z);
end
% Solve the saddle-point equations by damped Newton, starting from X:
%   g_r(u) = u_r*(Z_r + sum_k L_kr/(1-U_k)) - N_r = 0
u = X(:);
Zc = Z(:); Nc = N(:);
Uk = L*u;
if max(Uk) >= 1
    u = u * (1-1e-6)/max(Uk);
end
converged = false;
for it=1:200
    Uk = L*u;
    D = 1./(1-Uk);
    g = u.*(Zc + L'*D) - Nc;
    if norm(g) <= 1e-12*Ntot
        converged = true;
        break;
    end
    J = diag(Zc + L'*D) + (u*ones(1,R)).*(L'*(D.^2.*L));
    du = -J\g;
    alpha = 1;
    while any(u+alpha*du <= 0) || max(L*(u+alpha*du)) >= 1
        alpha = alpha/2;
        if alpha < 1e-12
            break;
        end
    end
    if alpha < 1e-12
        break;
    end
    u = u + alpha*du;
end
Uk = L*u;
D = 1./(1-Uk);
if converged && norm(u.*(Zc + L'*D) - Nc) <= 1e-8*Ntot
    us = u;   % exact saddle point
else
    us = X(:); % fallback: AQL/BS throughput (stationarity limits the damage)
end
% Assemble the expansion at us
Uk = L*us;
D = 1./max(GlobalConstants.FineTol, 1-Uk);
H = diag(Nc./us.^2) + L'*((D.^2).*L);
F = Zc'*us - sum(log(max(GlobalConstants.FineTol,1-Uk))) - Nc'*log(us);
lG = F - sum(log(us)) - (R/2)*log(2*pi) - 0.5*log(det(H)) + slcdemandfactor;
G = exp(lG);
end
