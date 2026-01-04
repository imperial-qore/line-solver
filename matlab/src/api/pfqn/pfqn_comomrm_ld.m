%{
%{
 % @file pfqn_comomrm_ld.m
 % @brief CoMoM for repairman model with load-dependent service rates.
%}
%}

%{
%{
 % @brief CoMoM for repairman model with load-dependent service rates.
 % @fn pfqn_comomrm_ld(L, N, Z, mu, options)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mu Load-dependent rate matrix (MxNt matrix).
 % @param options Solver options.
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
 % @return prob State probability distribution.
%}
%}
function [G,lG,prob] = pfqn_comomrm_ld(L,N,Z,mu, options)
% S: number of servers at the queueing stations
N = ceil(N);
atol = options.tol;
[M,R] = size(L);
Nt = sum(N);
Z = sum(Z,1);
if sum(Z) < GlobalConstants.Zero
    zset = false(M,1);
    for ist=1:M
        if norm(mu(ist,:)-(1:Nt),2) < atol
            zset(ist) = true;
        end
    end
    Z = L(zset,:);
    mu(zset,:)=[];
    L(zset,:)=[];
end

if sum(L) < GlobalConstants.Zero
    % model has only delays so it is trivial
    [G,lG] = pfqn_ca(L,N,Z);
    prob = zeros(sum(N)+1,1);    
    prob(end)=1;
    return
end
if nargin<4
    m=1;
end
[~,L,N,Z,lG0] = pfqn_nc_sanitize(zeros(1,R),L,N,Z,atol);
[M,R] = size(L);
if isempty(Z)
    if isempty(L)
        lG=lG0;
        G = exp(lG);
        prob = zeros(Nt+1,1); prob(1)=1;
        return
    end
    Z = zeros(1,R);
elseif isempty(L)
    L = zeros(1,R);
end
if M == 0
    G=exp(lG0);
    lG = lG0;
    prob = zeros(sum(N)+1,1);
    prob(end)=1;
    return
end

if M~=1
    line_error(mfilename,'The solver accepts at most a single queueing station.')
end

h = sparse(zeros(Nt+1,1)); h(Nt+1,1)=1;
scale = zeros(Nt,1);
nt = 0;
for r=1:R
    Tr = Z(r)*speye(Nt+1) + diag(sparse(L(r)*(Nt:-1:1)./mu(Nt:-1:1)),1);
    for nr=1:N(r)
        nt = nt + 1;
        h = Tr/nr * h;
        scale(nt) = abs(sum(sort(h))); % sort minimizes numerical issues
        h = abs(h)/scale(nt); % rescale so that |h|=1
    end
end

lG = lG0 + sum(log(scale));
G = exp(lG);
prob = h(end:-1:1)./G;
prob = prob/sum(prob);
end