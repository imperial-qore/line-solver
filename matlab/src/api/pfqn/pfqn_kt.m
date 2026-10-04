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
% A class with no jobs contributes a factor of 1 to G, but its saddle point is
% u_r -> 0, where N_r*log(u_r) and N_r/u_r^2 are indeterminate and lG comes back
% NaN. Solve the reduced model, as the self-looping branch below already does.
if any(N==0) && any(N>0)
    keep = N>0;
    [G,lG,X,Q] = pfqn_kt(L(:,keep), N(keep), Z(keep));
    return
end
[Morig,Rorig] = size(L);
% fix self-looping customers as they would yield Uk=1. Extracting u_r^N_r from
% 1/(1-U_ist) exactly leaves L(ist,r)^N_r and raises that station's factor to
% (1-V_ist)^-(1+N_r); classes looping at the SAME station share one factor
% (1+sum N) and contribute the multinomial (sum N)!/prod N_r!
slcdemandfactor = 0;
if Rorig>1
    isslc = false(1,Rorig);
    slcstation = zeros(1,Rorig);
    for r=1:Rorig
        % detect on the original rows: appended copies must not mask a later class
        if nnz(L(1:Morig,r))==1 && Z(r)==0
            isslc(r) = true;
            slcstation(r) = find(L(1:Morig,r)>0,1);
        end
    end
    for ist = unique(slcstation(isslc))
        grp = find(isslc & slcstation==ist);
        ntot = sum(N(grp));
        slcdemandfactor = slcdemandfactor + N(grp)*log(L(ist,grp))' ...
            + gammaln(ntot+1) - sum(gammaln(N(grp)+1));
        L = [L; repmat(L(ist,:),ntot,1)];
    end
    L(:,isslc)=[];
    Z(:,isslc)=[];
    N(:,isslc)=[];
end
[M,R] = size(L);
Ntot = sum(N);
if R==0 || Ntot==0 % nothing left to expand: the demand factors are exact
    lG = slcdemandfactor;
    G = exp(lG);
    X = zeros(1,R);
    Q = zeros(M,R);
    return;
end
% Solve the saddle-point equations by damped Newton:
%   g_r(u) = u_r*(Z_r + sum_k L_kr/(1-U_k)) - N_r = 0
% The start is closed form and pfqn_bs/pfqn_aql are consulted lazily, only if
% that fails. Seeding the iteration with an approximate MVA solve up front cost
% more than the Newton it seeded (96% of the runtime at M=8, R=4, N_r=32) and
% landed on the same saddle point to machine precision, so it is deferred to the
% failure path, where it still serves as a second start and as the last resort.
Zc = Z(:); Nc = N(:);
[u,converged] = pfqn_kt_newton(L,Zc,Nc,(0.5*N./max(sum(L,1),eps))',Ntot);
if converged
    us = u;   % exact saddle point
else
    if Ntot <= 4
        [X,Q] = pfqn_bs(L,N,Z); %#ok<ASGLU>
    else
        [X,Q] = pfqn_aql(L,N,Z); %#ok<ASGLU>
    end
    [u,converged] = pfqn_kt_newton(L,Zc,Nc,X(:),Ntot);
    if converged
        us = u;
    else
        us = X(:); % fallback: AQL/BS throughput (stationarity limits the damage)
    end
end
% X and Q now come from the saddle rather than from the approximate MVA solve:
% u* is the asymptotic MVA fixed point, so X_r = u*_r and Q_kr = u*_r L_kr/(1-U_k),
% which conserves the population exactly by the stationarity condition above.
X = us(:)';
Q = (L.*repmat(X,M,1))./repmat(max(GlobalConstants.FineTol,1-L*us),1,R);
% Assemble the expansion at us
Uk = L*us;
D = 1./max(GlobalConstants.FineTol, 1-Uk);
H = diag(Nc./us.^2) + L'*((D.^2).*L);
F = Zc'*us - sum(log(max(GlobalConstants.FineTol,1-Uk))) - Nc'*log(us);
% log|H| from the Cholesky factor: det(H) of an R x R Hessian leaves double range
% well before its logarithm does (it overflowed at R = 64, turning lG into -Inf)
[Hchol,cholfail] = chol(H);
if cholfail == 0
    logdetH = 2*sum(log(diag(Hchol)));
else % not numerically positive definite: fall back to the LU factors
    [~,Uh,~] = lu(H);
    logdetH = sum(log(abs(diag(Uh))));
end
lG = F - sum(log(us)) - (R/2)*log(2*pi) - 0.5*logdetH + slcdemandfactor;
G = exp(lG);
end

function [u,converged] = pfqn_kt_newton(L,Zc,Nc,u0,Ntot)
% Damped Newton on the saddle-point equations, from the start u0. Returns
% converged only if the residual also passes the outer 1e-8 check the caller
% used to apply, so a start that stalls is reported rather than silently used.
R = size(L,2);
u = u0;
Uk = L*u;
if max(Uk) >= 1
    u = u * (1-1e-6)/max(Uk);
end
converged = false;
for it=1:200 %#ok<NASGU>
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
D = 1./(1-L*u);
converged = converged && norm(u.*(Zc + L'*D) - Nc) <= 1e-8*Ntot;
end
