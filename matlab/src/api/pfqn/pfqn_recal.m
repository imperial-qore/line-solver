%{
%{
 % @file pfqn_recal.m
 % @brief RECAL (REcursive CALculation) method for normalizing constant.
%}
%}

%{
%{
 % @brief RECAL (REcursive CALculation) method for normalizing constant.
 % @fn pfqn_recal(L, N, Z, m0)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector (default: zeros).
 % @param m0 Initial multiplicity vector (default: ones).
 % @return G Normalizing constant.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [G,lG]=pfqn_recal(L,N,Z,m0)
% [G,logG]=PFQN_RECAL(L,N,Z,M0)
[M_original,R] = size(L);
if nargin<4
    m0=ones(1,M_original);
end
if nargin<3
    Z=zeros(1,R);
end

% Detect and consolidate replicated stations
[L, ~, ~, ~, mapping] = pfqn_unique(L);
[M,~] = size(L);

% Combine user-provided m0 with detected multiplicity
% For each unique station j, sum the m0 values of all original stations mapping to it
m0_combined = zeros(1, M);
for i = 1:M_original
    m0_combined(mapping(i)) = m0_combined(mapping(i)) + m0(i);
end
m0 = m0_combined;

Ntot = sum(N);
G_1 = ones(1,nchoosek(Ntot+(M+1)-1,Ntot));
G = G_1;
    I_1=multichoose(M+1,Ntot);
    n=0;
    for r=1:R
        for nr=1:N(r)
            n=n+1   ;
            I=multichoose(M+1,(Ntot+1)-(n+1));
            for i=1:size(I,1)
                m=I(i,:);
                mZ = m(1:M);
                G(i)=Z(r)*G_1(matchrow(I_1(:,1:M),mZ))/nr;
                for jst=1:M
                    m(jst)=m(jst)+1;
                    G(i)=G(i)+(m(jst)+m0(jst)-1)*L(jst,r)*G_1(matchrow(I_1,m))/nr;
                    m(jst)=m(jst)-1;
                end
            end
            I_1=I;
            G_1=G;
        end
    end
G = G(1);
lG = log(G);
end
