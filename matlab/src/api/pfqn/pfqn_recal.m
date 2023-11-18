function [G,lG]=pfqn_recal(L,N,Z,m0)
% [G,logG]=PFQN_RECAL(L,N,Z,M0)
[M,R] = size(L);
Ntot = sum(N);
G_1 = ones(1,nchoosek(Ntot+M-1,Ntot));
G = G_1;
if nargin<4
    m0=ones(1,M);
end
if nargin<3 || sum(Z)==0
    I_1=multichoose(M,Ntot);
    n=0;
    for r=1:R
        for nr=1:N(r)
            n=n+1   ;
            I=multichoose(M,(Ntot+1)-(n+1));
            for i=1:size(I,1)
                m=I(i,:);
                G(i)=(0);
                for j=1:M
                    m(j)=m(j)+1;
                    G(i)=G(i)+(m(j)+m0(j)-1)*L(j,r)*G_1(matchrow(I_1,m))/nr;
                    m(j)=m(j)-1;
                end
            end
            I_1=I;
            G_1=G;
        end
    end
else
    line_error('RECAL not yet available for models with think times.');
end
G = G(1);
lG = log(G);
end
