function met = spn_metrics(mdds, g, info)
% MET = SPN_METRICS(MDDS, G_L, INFO)
% Stationary measures of a product-form stochastic Petri net from the MDD-rec
% masses.
%
%   n(P_j) = sum_k k P(m_j = k)                    mean tokens
%   u(P_j) = 1 - P(m_j = 0)                        place utilization
%   u(T_j) = P(e_j >= 1)                           transition utilization
%   x(T_j) = sum_k min(k, c_j) W(T_j) P(e_j = k)   throughput
%   x(P_j) = sum_T I_j(T) x(T)                     tokens removed per unit time
%
% ONE DEVIATION FROM THE PAPER'S x(T_j), AND IT IS A GENERALISATION. The paper
% writes x(T_j) = sum_k k W(T_j) P(e_j = k), which is INFINITE-SERVER firing
% semantics -- every enabling set fires in parallel. LINE's own rate law is
% min(enabling degree, nmodeservers) * W(T), so c_j above is the mode's server
% count: c_j = 1 recovers single-server semantics, x = W(T) P(e >= 1), and
% c_j = Inf recovers the paper's formula exactly. Using the paper's form for a
% single-server mode would report a throughput that grows with the token
% population of a net whose transition can only fire one set at a time.
%
% The measures come out of ONE reachable set and ONE set of g_l, so they are
% mutually consistent by construction: no per-measure fixed point, no
% iteration.
%
% -- Input
% MDDS : MDD.toStruct of the reachable set built by SPN_MDD
% G_L  : 1 x K cell of per-level product-form factors g_l(v)
% INFO : the metadata SPN_MDD returned alongside the diagram
% -- Output
% MET : struct with fields G, tokens, placeUtil, placeTput, modeUtil,
%       modeTput and marginal, the per-level P(m_l = k)
%
% -- Reference
% S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
% product-form models of distributed systems with synchronisation", Future
% Generation Computer Systems 111 (2020) 475-490, Sec. 3.1 and Sec. 5.3.
%
% See also MDD_REC, MDD_REC_MARGINAL, SPN_REC_ENABLED, SPN_MDD.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

L = info.nplacelevels;
md = info.modes;

G = mdd_rec(mdds, g);
if ~(G > 0)
    line_error(mfilename, ['the normalising constant is not positive; the g_l passed do ' ...
        'not describe a product form over this reachable set']);
end

marginal = cell(1, L);
tokens = zeros(1, L);
placeUtil = zeros(1, L);
placeTput = zeros(1, L);
for l = 1:L
    pk = mdd_rec_marginal(mdds, g, l) / G;
    marginal{l} = pk;
    tokens(l) = (0:numel(pk) - 1) * pk(:);
    placeUtil(l) = 1 - pk(1);
end

E = numel(md);
modeUtil = zeros(1, E);
modeTput = zeros(1, E);
for e = 1:E
    en = spn_rec_enabled(mdds, g, md(e), L);
    modeUtil(e) = en.ge(2) / G;
    % W(T) is the scalar firing rate of the mode; a phase-type firing time has
    % no single rate, so its throughput is left to the phase-level marginal
    % rather than reported through this formula.
    if md(e).nph > 1
        line_error(mfilename, sprintf(['mode %d of node %d has a phase-type firing time, ' ...
            'whose throughput is not W(T) times an enabling probability; read it from the ' ...
            'phase-level marginal instead'], md(e).mode, md(e).trans));
    end
    % The formula below is W(T)*E[min(enabling degree, servers)], which is the
    % rate law only when no marking-dependent multiplier is in play. With one,
    % the firing rate is not a function of the enabling degree at all, so the
    % enabling-degree law is the wrong summary to take it from. The MARGINALS
    % above are unaffected -- they come from the product form, not the rates.
    if isfield(md, 'dep') && ~isempty(md(e).dep)
        line_error(mfilename, sprintf(['mode %d of node %d has a marking-dependent firing ' ...
            'rate, so its throughput is not W(T) times a function of the enabling degree and ' ...
            'cannot be read from the enabling-degree law. The token marginals are still ' ...
            'exact'], md(e).mode, md(e).trans));
    end
    rate = md(e).D1(1);
    k = (0:numel(en.eq) - 1);
    if isinf(md(e).srv), served = k; else, served = min(k, md(e).srv); end
    x = (served * rate) * en.eq(:) / G;
    modeTput(e) = x;
    for l = 1:L
        if md(e).enab(l) > 0
            placeTput(l) = placeTput(l) + md(e).enab(l) * x;
        end
    end
end

met = struct('G', G, 'tokens', tokens, 'placeUtil', placeUtil, 'placeTput', placeTput, ...
    'modeUtil', modeUtil, 'modeTput', modeTput, 'marginal', {marginal});
end
