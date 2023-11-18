function [lG] = pfqn_nrp(L,N,Z,alpha,options)
[M,R]=size(L);
Nt = sum(N); 
if sum(Z)>0
    L = [L;Z];
    alpha(end+1,1:Nt)=1:Nt;
end
if M==1 && sum(Z)==0
    [~,lG] = pfqn_gld(L,N,alpha);
    return
else
    Lmax = max(L); 
end
L = L./repmat(Lmax,size(L,1),1); % scale demands in [0,1]
x0 = zeros(1,R);
[~,~,lG] = laplaceapprox(@(x) infradius_hnorm(x, L, N,alpha),x0);
lG = lG + N*log(Lmax');
end