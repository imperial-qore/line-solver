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
if nargin<2
    c = zeros(1,M);
    mu = pfqn_fnc(alpha,c);
    if ~all(isfinite(mu)) % first retry with -1/2
        c = -0.5*ones(1,M);
        mu = pfqn_fnc(alpha,c);
    end
    dt = 0;
    it = 0;
    while ~all(isfinite(mu)) % randomize c if need be but unlikely
        it = it +1;
        dt = dt + 0.05;
        c = -1/2+dt;
        mu = pfqn_fnc(alpha,c);
        if c>=2
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