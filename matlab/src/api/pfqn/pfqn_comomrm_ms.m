%{
%{
 % @file pfqn_comomrm_ms.m
 % @brief CoMoM for multiserver repairman model.
%}
%}

%{
%{
 % @brief CoMoM for multiserver repairman model.
 % @fn pfqn_comomrm_ms(L, N, Z, m, S)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param m Replication factor (default: 1).
 % @param S Number of servers at queueing stations.
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
 % @return prob State probability distribution.
%}
%}
function [G,lG,prob] = pfqn_comomrm_ms(L,N,Z,m,S)
% m: replication factor
% S: number of servers at the queueing stations
[M,R] = size(L);
if M~=1
    line_error(mfilename,'The solver accepts at most a single queueing station.')
end
if nargin<4
    m=1;
end
atol = GlobalConstants.FineTol;
[~,L,N,Z,lG0] = pfqn_nc_sanitize(zeros(1,R),L,N,Z,atol);
Nt = sum(N);
if m>1
    mu = pfqn_mu_ms(Nt,m,S);
else
    mu = min(S,1:Nt);
end

h = zeros(Nt+1,1); h(Nt+1,1)=1;
scale = zeros(Nt,1);
nt = 0;
for r=1:R
    Tr = Z(r)*eye(Nt+1) + diag(L(r)*(Nt:-1:1)./mu(Nt:-1:1),1);
    for nr=1:N(r)
        nt = nt + 1;
        h = Tr/nr * h;
        scale(nt) = abs(sum(sort(h)));
        h = abs(h)/scale(nt); % rescale so that |h|=1
    end
end

lG = lG0 + sum(log(scale));
G = exp(lG);
prob = h(end:-1:1)./G;
prob = prob/sum(prob);
end