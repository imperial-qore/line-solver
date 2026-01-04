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

if nargin<4 || isempty(mu)
    mu = ones(m,sum(N)+1);
else
    mu = [1,mu(:)'];
end
if nargin<5
    m=1;
end
[~,R]=size(L);
% compute solution for [1,0,0,...,0]
p0 = zeros(sum(N)+1,1); p0(end)=1;
% compute the rest
tic;
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
if any(~isfinite(pk(1,1))) || ~isfinite(G)
    % todo
    lG = logsumexp(log(pk(1,:)));
elseif ~isfinite(G)
    lG = logsumexp(log(pk(1,:)));
else
    lG = log(G);
end
pk = pk/G;
pk= pk(end:-1:1);
%% test
Q=0;
% V=0;
for n=0:sum(N)
     psingle(n+1)= pk(n+1);
end
psingle=psingle/sum(psingle);
for n=1:sum(N)
     Q=Q + n * psingle(n+1);
%     V=V + (n^2-n) * psingle(n+1);  
end
% psingle
QN=double(Q)
[XNMVA,QNMVA,~,~,~,~,pik]=pfqn_mvald(repmat(L,m,1),N,Z,repmat(mu(:,2:end),m,1))
% VNMVA=0;
% for n=1:sum(N)
%     VNMVA=VNMVA+(n^2-n)*pik(1,:);
% end
% QNMVA = sum(QNMVA,2)
% QNTOTMVA=sum(QNMVA)
% [V,VNMVA]
end