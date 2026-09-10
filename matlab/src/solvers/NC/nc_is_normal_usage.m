function tf = nc_is_normal_usage(sn, form)
% TF = NC_IS_NORMAL_USAGE(SN, FORM)
%
% Is the closed model in NORMAL USAGE, the domain of the Mitra-McKenna PANACEA
% asymptotic expansion (J. ACM 33(3), 1986)?
%
% Normal usage asks that every queueing centre be able to absorb the load the
% think stations offer it: with rho_j0 = Ztot(j) the aggregate think demand of
% chain j, r_ij = L_ij / rho_j0 and mu_i(Ntot) the saturation rate,
%
%     alpha_i = 1 - (sum_j N_j r_ij) / mu_i(Ntot) > 0     for every centre i.
%
% Outside it the {phi(n)} series of the expansion DIVERGES, which is why
% PFQN_PANACEALD returns NaN there and PFQN_NCLD turns that NaN into a refusal
% rather than a numerical warning. This is a property of the demands and not of
% a construct the model declares, so it has no feature-registry name and cannot
% live in a feature set; NC_METHOD_REFUSAL is what carries it to the support
% gate.
%
% The rates are the ones SOLVER_NCLD would build: mu_i(n) = 1 for an ordinary
% single server, min(n,c) for a finite multiserver (the conversion runAnalyzer
% performs on the 'panald' arm), and the declared sn.lldscaling row when the
% model sets one. An infinite server is a think station and contributes to
% Ztot rather than to the centres.
%
% FORM selects which of the two shapes the demands take on the way to the
% expansion. 'lattice' (default) is the load-dependent one above, what
% PFQN_PANACEALD receives. 'seidmann' is the load-INDEPENDENT one of the
% 'pana' arm of PFQN_NC, which SOLVER_NC reaches through Seidmann's
% approximation: a c-server centre enters as demand L/c with mu_i(n) = 1 and its
% remainder L(c-1)/c is added to the think demand, so a multiserver model may be
% in normal usage in one form and not in the other, and the gate has to ask the
% form the run will take.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = true;
if nargin < 2 || isempty(form)
    form = 'lattice';
end
seidmann = strcmpi(form, 'seidmann');

% An open or mixed chain is refused earlier, by the closed-population feature
% set of the load-dependent evaluators; there is no rho_j0 to expand around.
if any(isinf(sn.njobs))
    return
end

[Lchain,~,~,~,Nchain] = sn_get_demands_chain(sn);
Nchain = Nchain(:)';
Nt = round(sum(Nchain(isfinite(Nchain))));
if Nt < 1
    return % the empty network: G = 1, nothing to expand
end

M = sn.nstations;
lld = sn.lldscaling;
if seidmann || isempty(lld)
    lld = ones(M, Nt);
    if ~seidmann
        for i = 1:M
            if isfinite(sn.nservers(i)) && sn.nservers(i) > 1
                lld(i,:) = min(1:Nt, sn.nservers(i));
            end
        end
    end
end
if size(lld,2) < Nt
    lld = [lld, repmat(lld(:,end), 1, Nt - size(lld,2))];
end

isIS = isinf(sn.nservers(:))';
Ztot = sum(Lchain(isIS,:), 1);
if isempty(Ztot)
    Ztot = zeros(1, numel(Nchain));
end
Lq = Lchain(~isIS, :);
if seidmann
    % Seidmann's split, as SOLVER_NC builds Lms and Zms before PFQN_NC
    cq = sn.nservers(~isIS);
    cq = cq(:);
    for i = 1:numel(cq)
        if isfinite(cq(i)) && cq(i) > 1
            Ztot = Ztot + Lq(i,:) * (cq(i) - 1) / cq(i);
            Lq(i,:) = Lq(i,:) / cq(i);
        end
    end
end
if any(Nchain > 0 & Ztot <= 0)
    % no think station on the route of a populated chain: the expansion
    % parameter rho_j0 is undefined
    tf = false;
    return
end

if isempty(Lq)
    return % no queueing centre: the expansion is the exact delay-only constant
end
muq = lld(~isIS, 1:Nt);
if any(muq(:) <= 0) || any(~isfinite(muq(:)))
    tf = false;
    return
end

r = zeros(size(Lq));
for j = 1:numel(Ztot)
    if Ztot(j) > 0
        r(:,j) = Lq(:,j) / Ztot(j);
    end
end
alpha = 1 - (r * Nchain(:)) ./ muq(:,Nt);
tf = min(alpha) > 0;
end
