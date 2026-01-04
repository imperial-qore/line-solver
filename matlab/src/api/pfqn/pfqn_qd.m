%{
%{
 % @file pfqn_qd.m
 % @brief Queue-Dependent (QD) approximate MVA solver.
%}
%}

%{
%{
 % @brief Queue-Dependent (QD) approximate MVA solver.
 % @fn pfqn_qd(L, N, ga, be, Q0)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param ga Gamma scaling function (default: ones).
 % @param be Beta scaling function (default: ones).
 % @param Q0 Initial queue length estimate (optional).
 % @return Q Mean queue lengths.
 % @return X System throughput.
 % @return U Utilization.
 % @return iter Number of iterations performed.
%}
%}
function [Q,X,U,iter] = pfqn_qd(L,N,ga,be,Q0)
[M,R]=size(L);

Q = zeros(M,R);
if nargin <3
    ga = @(A) ones(M,1);
end
if nargin <4
    be = @(A) ones(M,R);
end
if nargin < 5
    Q = L ./ repmat(sum(L,1),M,1) .* repmat(N,M,1);
else
    Q=Q0;
end
delta  = (sum(N) - 1) / sum(N);
deltar = (N - 1) ./ N;

Q_1 = Q*10;
tol = 1e-6;
iter = 0;
while max(max(abs(Q-Q_1))) > tol
    iter = iter + 1;
    Q_1 = Q;
    for k=1:M
        for r=1:R
            Ak{r}(k,1) = 1 + delta * sum(Q(k,:));
            Akr(k,r) = 1 + deltar(r) * Q(k,r);
        end
    end
    
    %    Q
    for r=1:R
        g = ga(Ak{r});
        b = be(Akr);
        for k=1:M
            C(k,r) = L(k,r) * g(k) * b(k,r) * (1 + delta * sum(Q(k,:)));
        end
        
        X(r) = N(r) / sum(C(:,r));
        
        for k=1:M
            Q(k,r) = X(r) * C(k,r);
            U(k,r) = L(k,r) * g(k) * b(k,r) * X(r);
        end
    end
end

end
