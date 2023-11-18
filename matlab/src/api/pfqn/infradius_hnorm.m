function y=infradius_hnorm(x,L,N,alpha)
M=size(L,1);
MU=0;
SIGMA=1;
Nt = sum(N);
beta = N/Nt;
t = normcdf(x,MU,SIGMA);
tb = sum(beta.*t,2);
h = @(x) pfqn_gld(sum(L.*repmat(exp(2*pi*1i*(t-tb)),M,1),2),Nt,alpha)* prod(normpdf(x,MU,SIGMA) ,2);

y = zeros(size(x,1),1);
for i=1:size(x,1)
    y(i) = real(h(x(i,:)));    
end
end
