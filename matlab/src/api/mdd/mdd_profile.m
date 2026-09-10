function mdd_profile(M, N)
% MDD_PROFILE(M, N)
% Profile the MDD reachability-set build for a single-class closed cyclic
% network with M stations and N jobs, printing the top hotspots by total time.
% Temporary developer harness (not part of the public API).
if nargin < 1, M = 6; end
if nargin < 2, N = 20; end
domain = (N + 1) * ones(1, M);
init = zeros(1, M); init(1) = N;

profile clear
profile on
t0 = tic;
mdd = mdd_reachset(domain, init, @(s) i_cyc(s));
elapsed = toc(t0);
profile off

s = mdd.stats();
fprintf('M=%d N=%d  |S|=%d  nodes=%d  build=%.2fs\n', M, N, s.numStates, s.numNodes, elapsed);

info = profile('info');
ft = info.FunctionTable;
[~, ix] = sort([ft.TotalTime], 'descend');
fprintf('%10s %10s %10s   %s\n', 'total(s)', 'self(s)', 'calls', 'function');
for i = 1:min(18, numel(ix))
    f = ft(ix(i));
    self = f.TotalTime - sum([f.Children.TotalTime]);
    fprintf('%10.3f %10.3f %10d   %s\n', f.TotalTime, self, f.NumCalls, f.FunctionName);
end
end

function T = i_cyc(s)
M = numel(s);
T = zeros(M, M);
m = 0;
for i = 1:M
    if s(i) > 0
        j = mod(i, M) + 1;
        t = s; t(i) = t(i) - 1; t(j) = t(j) + 1;
        m = m + 1;
        T(m, :) = t;
    end
end
T = T(1:m, :);
end
