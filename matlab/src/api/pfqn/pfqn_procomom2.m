%{
%{
 % @file pfqn_procomom2.m
 % @brief Product-form CoMoM for 2-station repairman model (queue + delay).
%}
%}

%{
%{
 % @brief Product-form CoMoM for 2-station repairman model (queue + delay).
 % @fn pfqn_procomom2(L, N, Z, mu, m)
 % @param L Service demand vector.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mu Load-dependent rates (optional).
 % @param m Replication factor (default: 1).
 % @return pk Marginal state probabilities.
 % @return lG Logarithm of normalizing constant.
 % @return G Normalizing constant.
 % @return T Transfer matrices.
 % @return F Product transfer matrix.
 % @return B Combined transfer matrix.
%}
%}
function [pk,lG,G,T,F,B]=pfqn_procomom2(L,N,Z,mu,m)
% Marginal state probabilities for the queue in a model consisting of a
% queueing station and a delay station only.

% m must be defaulted BEFORE the mu block, which reads it. With the two blocks
% in the opposite order every call with fewer than five arguments died on an
% undefined m ("Not enough input arguments"), so the 3- and 4-argument forms
% advertised in the header above were dead.
if nargin<5
    m=1;
end
if nargin<4 || isempty(mu)
    mu = ones(m,sum(N)+1);
else
    mu = [1,mu(:)'];
end
[~,R]=size(L);
% compute solution for [1,0,0,...,0]
p0 = zeros(sum(N)+1,1); p0(end)=1;
% compute the rest
for r=1:R
    % generate F2r matrix
    T{r} = sparse(1+sum(N),1+sum(N));
    for n=sum(N):-1:1
        row = sum(N)-n+1;
        T{r}(row,row) = Z(r);
        T{r}(row,row+1) = (n+m-1)*L(r)/mu(1+n);
    end
    T{r}(sum(N)+1,sum(N)+1) = Z(r);    
end
F = eye(sum(N)+1);
B = eye(sum(N)+1);
for r=1:R
    F = F*T{r}^N(r)/factorial(N(r));
    B = B*T{r};
end
pk = (F*p0)';
G = sum(pk);
% The former middle branch (elseif ~isfinite(G)) was unreachable: the first
% condition already contains ~isfinite(G).
if any(~isfinite(pk(1,1))) || ~isfinite(G)
    lG = logsumexp(log(pk(1,:)));
else
    lG = log(G);
end
pk = pk/G;
pk= pk(end:-1:1);
end
