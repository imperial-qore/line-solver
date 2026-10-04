clear node jobclass

% Open Whittle network: a bandwidth-sharing model in which one route holds
% SEVERAL links at once, which no per-station rate scaling can express. This is
% the 2-link linear network: route 1 crosses both links, routes 2 and 3 use one
% link each.
%
%   A = [1 1 0;      link 1 carries routes 1 and 2
%        1 0 1]      link 2 carries routes 1 and 3
%
% Capacity is shared by BALANCED FAIRNESS, whose rates come from the recursion
%   Phi(n) = max_l (1/C_l) sum_{s in l} Phi(n-e_s),   x_s(n) = Phi(n-e_s)/Phi(n).
% That construction satisfies the Whittle balance property by design, so the
% stationary law is pi(n) ~ Phi(n) prod rho_s^{n_s} and is insensitive.
%
% Each route is modelled as its own PS queue fed by its own open class, so the
% queue populations ARE the coordinates of the Whittle state n.
A = [1 1 0; 1 0 1];
C = [1; 1];
nu = [0.20 0.30 0.30];  % flow arrival rates
mu = [1.00 1.00 1.00];  % 1/mean file size
cutoff = 6;             % per-route truncation of the open state space
%%
[Phi, base] = bfBalanceFunction(A, C, cutoff*ones(3,1));

bwmodel = Network('model');
node{1} = Source(bwmodel, 'Source');
node{2} = Queue(bwmodel, 'Route1', SchedStrategy.PS);
node{3} = Queue(bwmodel, 'Route2', SchedStrategy.PS);
node{4} = Queue(bwmodel, 'Route3', SchedStrategy.PS);
node{5} = Sink(bwmodel, 'Sink');
for s = 1:3
    jobclass{s} = OpenClass(bwmodel, sprintf('Route%dFlows', s));
    node{1}.setArrival(jobclass{s}, Exp(nu(s)));
end
for s = 1:3
    for t = 1:3
        if s == t
            node{1+t}.setService(jobclass{s}, Exp(mu(s)));
        else
            node{1+t}.setService(jobclass{s}, Disabled.getInstance());
        end
    end
end
P = bwmodel.initRoutingMatrix();
for s = 1:3
    P{jobclass{s},jobclass{s}}(node{1},node{1+s}) = 1;
    P{jobclass{s},jobclass{s}}(node{1+s},node{5}) = 1;
end
bwmodel.link(P);

sn = bwmodel.getStruct();
idx = zeros(1,3);
for s = 1:3
    idx(s) = sn.nodeToStation(bwmodel.getNodeIndex(sprintf('Route%d', s)));
end
M = sn.nstations; K = sn.nclasses;
% The third argument is the per-slot open-class truncation used when phi is
% materialized onto the JSON wire; it matches options.cutoff below, which is also
% the range over which the balance function Phi was built.
bwmodel.setGlobalDependence(@(n) bfSpread(n, idx, Phi, base, M, K), 1.0, cutoff);

options = SolverCTMC.defaultOptions;
options.cutoff = cutoff;
bwAvgTableCTMC = CTMC(bwmodel, options).getAvgTable

% Cross-check against the closed-form product form. Truncating a reversible
% chain preserves the conditional law, so this agrees to solver precision
% rather than to a truncation-limited tolerance.
En = bfMeans(Phi, base, nu./mu);
line_printf('\nproduct form E[n] = [%.6f %.6f %.6f]\n', En);

model = bwmodel; % for test compatibility

function v = bfSpread(n, idx, Phi, base, M, K)
S = numel(idx);
np = zeros(S,1);
for s = 1:S
    np(s) = n(idx(s), s);
end
v = ones(M,K);
if all(np == 0)
    return
end
den = Phi(bfIndex(np, base));
for s = 1:S
    if np(s) > 0
        m = np; m(s) = m(s)-1;
        v(idx(s), s) = Phi(bfIndex(m, base)) / den;
    else
        v(idx(s), s) = 0;
    end
end
end

function [Phi, base] = bfBalanceFunction(A, C, cutoff)
S = size(A,2);
base = cutoff(:) + 1;
ns = prod(base);
states = bfStates(base, ns, S);
[~,ord] = sort(sum(states,2));
Phi = zeros(ns,1);
Phi(bfIndex(zeros(S,1), base)) = 1;
for kk = 1:ns
    k = ord(kk);
    n = states(k,:)';
    if all(n == 0)
        continue
    end
    best = 0;
    for l = 1:size(A,1)
        acc = 0;
        for s = 1:S
            if A(l,s) > 0 && n(s) > 0
                m = n; m(s) = m(s)-1;
                acc = acc + Phi(bfIndex(m, base));
            end
        end
        best = max(best, acc/C(l));
    end
    Phi(k) = best;
end
end

function En = bfMeans(Phi, base, rho)
S = numel(base);
ns = prod(base);
states = bfStates(base, ns, S);
w = zeros(ns,1);
for k = 1:ns
    w(k) = Phi(k) * prod(rho(:)'.^states(k,:));
end
w = w / sum(w);
En = w' * states;
end

function states = bfStates(base, ns, S)
states = zeros(ns,S);
for k = 1:ns
    rem = k-1;
    for s = 1:S
        states(k,s) = mod(rem, base(s));
        rem = floor(rem/base(s));
    end
end
end

function i = bfIndex(n, base)
i = 1; mult = 1;
for s = 1:numel(n)
    i = i + n(s)*mult;
    mult = mult*base(s);
end
end
