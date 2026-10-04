function [gamma, lambda_cache, Rcost] = da_cache_isolate(ch, lambda)
% [GAMMA, LAMBDA_CACHE, RCOST] = DA_CACHE_ISOLATE(CH, LAMBDA)
%
% Isolated-cache input construction for DA methods: builds the per-class,
% per-item, per-list arrival rates LAMBDA_CACHE and the access cost
% matrices RCOST from the cache node parameters CH (sn.nodeparam of the
% cache node) and the current class arrival rates LAMBDA (1 x nclasses),
% then computes the access factors GAMMA via cache_gamma_lp.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

m = ch.itemcap;
n = ch.nitems;
h = length(m);
u = length(lambda);
lambda_cache = zeros(u, n, h);

for v = 1:u
    for k = 1:n
        for l = 1:(h+1)
            if ~isnan(ch.pread{v})
                lambda_cache(v, k, l) = lambda(v) * ch.pread{v}(k);
            end
        end
    end
end

Rcost = ch.accost;
if isempty(Rcost)
    % Default linear cache routing: items flow from list l to list l+1
    Rcost = cell(u, n);
    for v = 1:u
        for k = 1:n
            Rmat = diag(ones(1, h), 1);
            Rmat(h+1, h+1) = 1;
            Rcost{v, k} = Rmat;
        end
    end
end
gamma = cache_gamma_lp(lambda_cache, Rcost);
end
