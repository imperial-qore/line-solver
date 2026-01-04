function ci = mtrace_bootstrap(T,A,resamples)

if nargin == 2
    resamples = 1000;
end

% number of samples
N = length(T);
% target block length
BL = 50;
% number of blocks
BN = floor(N / BL);
R = mod(N, BN);
% actual block lengths
L = ones(BN,1) * floor(N / BN);
L(1:R) = L(1:R) + 1;
% block indices
I = (1:BN)';
% compute confidence intervals
statfunc = @(idx) bf(idx);
ci = bootci(resamples, statfunc, I);


function STATS = bf(idx)
    [t,a] = blocks(T,A,idx,L);
    p = mtrace_pc(t,a);
    B = mtrace_moment(t,a,1,0,1);
    F = mtrace_moment(t,a,1,1,1);
    S = mtrace_sigma(t, a);
    STATS = [p', B', F', S(:)'];
end

function [t,a] = blocks(T,A,idx,L)
    n = sum(L(idx));
    t = zeros(n,1);
    a = zeros(n,1);
    last = cumsum(L);
    first = last - L + 1;
    next = 1;
    for i = 1:length(idx)
        j = idx(i);
        d1 = next;
        d2 = next + L(j) - 1;
        s1 = first(j);
        s2 = last(j);
        t(d1:d2) = T(s1:s2);
        a(d1:d2) = A(s1:s2);
        next = d2 + 1;
    end
end

end