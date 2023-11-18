function prob=cache_prob_asy(gamma,m)
% FPI method
[n,h]=size(gamma);
xi = cache_xi_fp(gamma,m);
prob = zeros(n,h);
for i=1:n
    prob(i,1) = 1/(1+gamma(i,:)*xi(:));
    prob(i,2:(1+h)) = gamma(i,:)*xi(:) ./ (1+gamma(i,:)*xi(:));
end
end
