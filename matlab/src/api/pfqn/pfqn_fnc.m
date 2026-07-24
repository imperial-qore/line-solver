%{
%{
 % @file pfqn_fnc.m
 % @brief Generate load-dependent rates for functional server model f(n)=n+c.
%}
%}

%{
%{
 % @brief Generate load-dependent rates for functional server model f(n)=n+c.
 % @fn pfqn_fnc(alpha, c)
 % @param alpha Rate parameters (Mx N matrix).
 % @param c Constant offset parameter (default: auto-determined).
 % @return mu Load-dependent service rates.
 % @return c Determined offset constant.
%}
%}
function [mu,c] = pfqn_fnc(alpha,c)
% generate rates for functional server f(n)=n+c
M = size(alpha,1);
if size(alpha,2) == 0
    % No rate columns: the caller shifted a single-column mu with
    % pfqn_mushift, which returns M x (N-1) and so yields zero columns when the
    % total population is 1. There is no functional-server rate to build, and
    % the population-N-1 subproblems the caller then forms are empty (G = 1).
    % Without this guard alpha(ist,1) below indexes an empty array and
    % SolverNC(model,'method','exact') fails on EVERY closed model whose total
    % population is 1.
    mu = zeros(M,0);
    c = zeros(1,M);
    return
end
if nargin<2
    c = zeros(1,M);
    mu = pfqn_fnc(alpha,c);
    % all(isfinite(mu)) reduces per COLUMN for M>1, so the ladder fired only
    % when EVERY population column held a non-finite rate; the intent is any.
    if any(~isfinite(mu(:))) % first retry with -1/2
        c = -0.5*ones(1,M);
        mu = pfqn_fnc(alpha,c);
    end
    dt = 0;
    it = 0;
    while any(~isfinite(mu(:))) % randomize c if need be but unlikely
        it = it +1;
        dt = dt + 0.05;
        % c must stay a 1xM vector: the recursive call indexes c(ist) per
        % station, so assigning a scalar here made this retry path fail with
        % "Index exceeds array bounds" for any M > 1.
        c = (-1/2+dt)*ones(1,M);
        mu = pfqn_fnc(alpha,c);
        if (-1/2+dt) >= 2
            break
        end
    end
    return
end
N = length(alpha(1,:));
mu = zeros(M,N);
for ist=1:M
    mu(ist,1) = alpha(ist,1)/(1+c(ist));
    alphanum = sparse(zeros(N,N));
    alphaden = sparse(zeros(N,N));
    for n=2:N
        alphanum(n,1) = alpha(ist,n);
        alphaden(n,1) = alpha(ist,n-1);
        for k=2:(n-1)
            alphanum(n,k) = alphanum(n,k-1) * alpha(ist,n-k+1);
            alphaden(n,k) = alphaden(n,k-1) * alpha(ist,n-k);
        end
    end
    for n=2:N
        rho = 0;
        muden = 1;
        for k=1:(n-1)
            muden = muden * mu(ist,k);
            rho = rho+(alphanum(n,k)-alphaden(n,k)) / muden;
        end
        mu(ist,n) = alphanum(n,n-1)*alpha(ist,1)/muden;
        mu(ist,n) = mu(ist,n)/(1-rho);
    end
end
mu(isnan(mu)) = Inf;
mu(abs(mu)>1e15) = Inf;
for ist=1:M
    if any(isinf(mu(ist,:)))
        s = min(find(isinf(mu(ist,:))));
        mu(ist,s:end)=Inf;
    end
end
%mu(mu==0) = Inf;
end