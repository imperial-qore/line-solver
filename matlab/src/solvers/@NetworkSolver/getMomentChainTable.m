function [MomentChainTable, mom] = getMomentChainTable(self, order)
% GETMOMENTCHAINTABLE Exact higher moments of the per-chain queue length.
%
% [MOMENTCHAINTABLE, MOM] = GETMOMENTCHAINTABLE(SELF) returns a table with one
% row per (Station, Chain) giving the moments of the queue length of that
% chain at that station, Q_(i,c) = sum_(r in chain c) n(i,r):
%   QLen, QLenVar, QLenSCV
%
% This is the chain-level analogue of getAvgChainTable, and it sits between the
% two other moment tables: getMomentTable is per class, getMomentStationTable
% is per station total, and this one is per chain, i.e. per group of classes
% that circulate together.
%
% [..] = GETMOMENTCHAINTABLE(SELF, ORDER) selects which moment orders to
% report. ORDER is a set: a scalar k is read as 1:k, "everything up to order
% k"; an explicit vector selects exactly those orders.
%   1        the mean only:             QLen
%   2        (default) mean and second moment: QLen, QLenVar, QLenSCV
%   3        also adds QLenM3 and QLenSkew
%
% Unlike the per-class table, order 3 IS available here. All three tables are
% the same recursion under different groupings of the classes: the generating
% parameter scales the service times of a class subset T at a station, and the
% moments it produces are those of sum_(r in T) n(i,r). T = {r} gives
% getMomentTable, T = chain gives this table, T = all classes gives
% getMomentStationTable. That is Theorem 1 of Akyildiz and Strelen; Strelen's
% own x_i is the last case.
%
% The ALGORITHM is chosen by the solver's method, set at construction, not by
% an argument here; see getMomentStationTable. A Linearizer-family method
% approximates the per-station totals only, so it cannot express a per-chain
% grouping and is rejected here unless every class already sits in one chain, in
% which case the chain IS the station total.
%
% Restricted to closed, single-server models, which is the scope of
% pfqn_sens_mom.
%
% MOM is the underlying pfqn_sens_mom struct. Its .Cov is (M x C x M x C) and
% carries the cross-chain and cross-station covariances this table does not
% show.
%
% Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
% Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
% Communications 39(6):828-832, 1991, Theorem 1; J. C. Strelen, "Moment
% Analysis for Closed Queuing Networks and its Linearizer", Performance
% Evaluation 11:127-142, 1990, equation (3.2).
%
% See also: getMomentTable, getMomentStationTable, getAvgChainTable.

if nargin < 2 || isempty(order)
    order = 2;
end
order = validateMomentOrder(order, 3);
% see _kb/06-solver-catalog.md for rationale (method fixed at construction)
method = self.getOptions.method;

sn = self.model.getStruct();
R = sn.nclasses;
N = sn.njobs;

[~, D, Np, Z, ~, Ssrv, ~] = sn_get_product_form_params(sn);
queueIndices = find(sn.nodetype == NodeType.Queue);
Mq = numel(queueIndices);
Ztot = sum(Z, 1);

% see _kb/06-solver-catalog.md for rationale (product-form identity Cov=L dQ/dL)
if ~sn_has_product_form(sn)
    line_error(mfilename, 'getMomentChainTable requires a product-form model: the moment identity Cov[n,n] = L dQ/dL holds only under product form, so no correct value exists here. Use SolverCTMC (exact distribution) or SolverLDES with setReward for the moments of a non-product-form model.');
end
if any(isinf(N))
    line_error(mfilename, 'getMomentChainTable supports closed models only. Per-class second moments of an open or mixed model are available from getMomentTable.');
end
if any(Ssrv > 1)
    line_error(mfilename, 'getMomentChainTable supports single-server stations only: the higher-moment recursion of the reference is stated for load-independent stations.');
end

% the chain of each class; sn.chains is (nchains x nclasses)
groups = zeros(1, R);
for r = 1:R
    c = find(sn.chains(:, r), 1);
    if isempty(c)
        line_error(mfilename, sprintf('class %s belongs to no chain', sn.classnames{r}));
    end
    groups(r) = c;
end
% pfqn_sens_mom requires the group labels to be consecutive from 1, so drop any
% chain that holds no class rather than leaving a hole in the numbering
[used, ~, compact] = unique(groups);
groups = compact(:)';
Cg = numel(used);

% see _kb/06-solver-catalog.md for rationale (Linearizer cannot express per-chain grouping)
if isLinearizerMethod(method)
    if Cg > 1
        line_error(mfilename, sprintf('the solver method ''%s'' approximates the per-station totals, so it cannot produce a per-chain grouping of %d chains. Use an exact method, or getMomentStationTable.', method, Cg));
    end
    mom = pfqn_sens_linearizer(D, Np, Ztot);
elseif isExactMvaMethod(method)
    mom = pfqn_sens_mom(D, Np, Ztot, ones(1,Mq), groups);
else
    % Finite-difference oracle over one chain's classes (Akyildiz-Strelen Thm 1
    % class subset T). see _kb/06-solver-catalog.md for rationale
    mom = chainMomentsByFiniteDifference(self, sn, queueIndices, groups, Cg);
end

Station = {}; Chain = {};
QLen = []; QLenVar = []; QLenSCV = []; QLenM3 = []; QLenSkew = [];
for ist = 1:Mq
    for g = 1:Cg
        classesOf = find(groups == g);
        if all(all(D(ist, classesOf) <= 0))
            continue;   % no class of this chain visits this station
        end
        Station{end+1, 1} = sn.nodenames{queueIndices(ist)}; %#ok<AGROW>
        Chain{end+1, 1} = sprintf('Chain%d', used(g));       %#ok<AGROW>
        QLen(end+1, 1) = mom.m(ist, g);          %#ok<AGROW>
        QLenVar(end+1, 1) = mom.Var(ist, g);     %#ok<AGROW>
        if mom.m(ist, g) > 0
            QLenSCV(end+1, 1) = mom.Var(ist, g) / mom.m(ist, g)^2; %#ok<AGROW>
        else
            QLenSCV(end+1, 1) = NaN;             %#ok<AGROW>
        end
        QLenM3(end+1, 1) = mom.M3(ist, g);       %#ok<AGROW>
        QLenSkew(end+1, 1) = mom.Skew(ist, g);   %#ok<AGROW>
    end
end

vars = {Station, Chain};
names = {'Station', 'Chain'};
if any(order == 1)
    vars{end+1} = QLen;     names{end+1} = 'QLen';
end
if any(order == 2)
    vars{end+1} = QLenVar;  names{end+1} = 'QLenVar';
    vars{end+1} = QLenSCV;  names{end+1} = 'QLenSCV';
end
if any(order == 3)
    vars{end+1} = QLenM3;   names{end+1} = 'QLenM3';
    vars{end+1} = QLenSkew; names{end+1} = 'QLenSkew';
end
MomentChainTable = table(vars{:}, 'VariableNames', names);
end

% =========================================================================
function tf = isExactMvaMethod(method)
% Methods whose means are the exact MVA recursion; see getMomentStationTable.
tf = any(strcmpi(method, {'default', 'mva', 'exact'}));
end

% =========================================================================
function mom = chainMomentsByFiniteDifference(self, sn, queueIndices, groups, Cg)
% Per-chain moments from ANY solver and method, by central differences of that
% method's OWN mean queue lengths; see solveMeansForStruct.
%
% The parameter y_(i,g) scales the demands of group g's classes at station i,
% which is the class-subset parameter T of Akyildiz-Strelen Theorem 1; the
% moments it generates are those of Q_(i,g) = sum_(r in g) n(i,r). Since
% D(i,r) = visits(i,r)/rate(i,r), the perturbation is applied to the rates of
% that group's classes at that station alone, leaving the other classes' demands
% at the same station untouched. That per-class granularity is what separates
% this from getMomentStationTable, which scales a whole column.
M = numel(queueIndices);
h = 1e-4;
qst = sn.nodeToStation(queueIndices);
P = M*Cg;
pidx = zeros(M,Cg);
p = 0;
for i = 1:M
    for g = 1:Cg
        p = p + 1; pidx(i,g) = p;
    end
end
m0 = solveGroupTotals(self, sn, qst, groups, Cg);
dm = zeros(M,Cg,M,Cg); d2m = zeros(M,Cg);
for hi = 1:M
    for hg = 1:Cg
        snp = scaleGroupDemands(sn, qst(hi), groups, hg, 1+h);
        snm = scaleGroupDemands(sn, qst(hi), groups, hg, 1-h);
        mp = solveGroupTotals(self, snp, qst, groups, Cg);
        mm = solveGroupTotals(self, snm, qst, groups, Cg);
        for i = 1:M
            for g = 1:Cg
                dm(i,g,hi,hg) = (mp(i,g) - mm(i,g)) / (2*h);
            end
        end
        d2m(hi,hg) = (mp(hi,hg) - 2*m0(hi,hg) + mm(hi,hg)) / h^2;
    end
end
mom = packChainFiniteDifference(m0, dm, d2m, Cg);
end

% =========================================================================
function sn2 = scaleGroupDemands(sn, station, groups, g, factor)
% Scale the demands of group g's classes at one station by FACTOR, via rates.
sn2 = sn;
cls = find(groups == g);
sn2.rates(station, cls) = sn.rates(station, cls) / factor;
end

% =========================================================================
function mg = solveGroupTotals(self, sn, qst, groups, Cg)
% Per-(station,group) mean queue lengths under the solver's own method.
QN = solveMeansForStruct(self, sn);
M = numel(qst);
mg = zeros(M,Cg);
for i = 1:M
    for g = 1:Cg
        mg(i,g) = sum(QN(qst(i), groups == g));
    end
end
end

% =========================================================================
function mom = packChainFiniteDifference(m, dm, d2m, Cg)
% (3.2), applied to numerically obtained derivatives, per (station,group).
M = size(m,1);
mom.m = m; mom.d2m = d2m;
flat = reshape(dm, M*Cg, M*Cg);
mom.CovAsym = max(max(abs(flat - flat.')));
flat = (flat + flat.')/2;
if Cg == 1
    mom.Cov = reshape(flat, M, M); mom.dm = reshape(dm, M, M);
else
    mom.Cov = reshape(flat, M, Cg, M, Cg); mom.dm = dm;
end
Var = zeros(M,Cg); M2 = zeros(M,Cg); M3 = zeros(M,Cg); Skew = zeros(M,Cg);
for i = 1:M
    for g = 1:Cg
        d1 = dm(i,g,i,g);
        Var(i,g) = d1;
        M2(i,g) = d1 + m(i,g)^2;
        M3(i,g) = d2m(i,g) + (1 + 3*m(i,g))*d1 + m(i,g)^3;
        mu3 = M3(i,g) - 3*m(i,g)*M2(i,g) + 2*m(i,g)^3;
        if Var(i,g) > 0
            Skew(i,g) = mu3 / Var(i,g)^1.5;
        else
            Skew(i,g) = NaN;
        end
    end
end
mom.Var = Var; mom.M2 = M2; mom.M3 = M3; mom.Skew = Skew;
end

% =========================================================================
function tf = isLinearizerMethod(method)
% True for the Linearizer family of solver methods; see getMomentStationTable.
tf = any(strcmpi(method, {'lin', 'amva.lin', 'egflin', 'gflin'}));
end

% =========================================================================
function order = validateMomentOrder(order, maxorder)
% ORDER is a set of moment orders. A scalar k is shorthand for 1:k, so that
% getMomentChainTable(2) means "up to the second moment" and not "the second
% moment alone"; a vector of two or more entries is taken literally.
%
% Consequence of MATLAB's isscalar: a one-element vector IS a scalar, so [2]
% takes the 1:k path and yields [1 2]. "The second moment alone" is therefore
% not expressible, which is deliberate: a variance with no mean beside it is not
% a useful table, and [2 3] remains available for the higher orders.
%
% Non-integers are rejected on BOTH paths. Rounding them silently would accept
% [1 2.5] as [1 3], i.e. answer a question that was not asked.
if ~isnumeric(order) || isempty(order) || any(~isfinite(order(:)))
    line_error(mfilename, sprintf('order must be an integer in 1..%d, or a vector of such integers.', maxorder));
end
if any(order(:) ~= round(order(:))) || any(order(:) < 1) || any(order(:) > maxorder)
    line_error(mfilename, sprintf('order must be an integer in 1..%d, or a vector of such integers.', maxorder));
end
if isscalar(order)
    order = 1:order;
    return;
end
order = unique(order(:)');
end
