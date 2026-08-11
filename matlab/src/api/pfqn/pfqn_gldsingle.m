%{
%{
 % @file pfqn_gldsingle.m
 % @brief Exact normalizing constant for single-class load-dependent models.
%}
%}

%{
%{
 % @brief Exact normalizing constant for single-class load-dependent models.
 % @fn pfqn_gldsingle(L, N, mu, options)
 % @param L Service demand vector (Mx1).
 % @param N Population (scalar).
 % @param mu Load-dependent rate matrix (MxN).
 % @param options Solver options.
 % @return lG Logarithm of normalizing constant.
 % @return G Normalizing constant.
%}
%}
function [lG,G]=pfqn_gldsingle(L,N,mu,options)
% G=PFQN_GLDSINGLE(L,N,MU)

if nargin<4
    options = [];
end

[M,R]=size(L);
if R>1
    line_error(mfilename,'multiclass model detected. pfqn_gldsingle is for single class models.');
end
Nscal = N(1); % codegen: ensure scalar loop bound

% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
useLog = isreal(L) && isreal(mu) && all(L(:)>=0) && all(mu(:)>0);

if useLog
    % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
    lg = -Inf(M+1, Nscal+1, Nscal+2);
    % lg(0+1,n+1,1+1) stays -Inf for n>=1: no station can hold n>=1 jobs.
    lL = log(L);       % -Inf where the demand is zero
    lmu = log(mu);     % +Inf where the rate is infinite, zeroing the term
    for m=1:M
        for tm=1:(Nscal+1)
            lg(m +1,0 +1,tm +1)=0; % log(1): zero jobs
        end
        for n=1:Nscal
            for tm=1:(Nscal-n+1)
                a = lg(m-1 +1, n +1, 1 +1);
                b = lL(m) + lg(m +1, n-1 +1, tm+1 +1) - lmu(m,tm);
                % pairwise log-sum-exp of a and b, stable when either is -Inf
                if a > b
                    if b == -Inf
                        lg(m +1, n +1, tm +1) = a;
                    else
                        lg(m +1, n +1, tm +1) = a + log1p(exp(b-a));
                    end
                else
                    if a == -Inf
                        lg(m +1, n +1, tm +1) = b;
                    else
                        lg(m +1, n +1, tm +1) = b + log1p(exp(a-b));
                    end
                end
            end
        end
    end
    lG = lg(M +1,Nscal +1,1 +1);
    G = exp(lG);
else
    g = zeros(M+1, Nscal+1, Nscal+2);
    for n=1:Nscal
        g(0 +1,n +1, 1 +1)=0;
    end
    for m=1:M
        for tm=1:(Nscal+1)
            g(m +1,0 +1,tm +1)=1;
        end
        for n=1:Nscal
            for tm=1:(Nscal-n+1)
                g(m +1, n +1, tm +1)= g(m-1 +1, n +1, 1 +1)+L(m)*g(m +1, n-1 +1, tm+1 +1)/mu(m,tm);
            end
        end
    end
    G = g(M +1,Nscal +1,1 +1);
    lG = log(G);
end
end