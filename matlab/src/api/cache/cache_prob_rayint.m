function prob = cache_prob_rayint(gamma,m, lE)
[n,h]=size(gamma);
if nargin < 3
[~, lE] = cache_rayint(gamma,m);
end
for i=1:n
    for j=1:h
        [~, lEi] = cache_rayint(gamma(setdiff(1:n,i),:),oner(m,j));
        prob(i,1+j) = m(j) * gamma(i,j) * exp(lEi-lE);
    end
    prob(i,1) = abs(1 - sum(prob(i,2:end)));
end
end