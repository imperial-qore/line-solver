function poc_ansim_prio()
% POC_ANSIM_PRIO  Priority-station comparison for the ansim surrogate.
%
% Model: Delay + Queue1 under PSPRIO (only the most urgent non-empty priority
% group is served, PS inside the group), D=(2,3), Z=(4,6), priorities (0,1),
% N1=N2=N, service law at the queue swept over exp / hyper-exp (SCV=4) /
% hypo-exp (SCV=0.5) at fixed mean. Reference: exact PSPRIO CTMC. Error is L1
% over the whole station-by-class queue-length matrix.
%
% Competitors fall into three groups.
%   (a) DPS-emulating: the discipline is replaced by DPS weights
%       w_r = W^(maxprio-prio_r). SolverFLD and SolverNC can take no other
%       shape. Floor: dps-gap = exact DPS-W CTMC vs exact PSPRIO CTMC.
%   (b) HOL-relaxing: non-preemptive priority. SolverMVA's HOL arm with both
%       np_priority settings ('cl' = Chandy-Lakshmi in the Eager-Lipscomb
%       arrival-instant form, 'shadow' = Sevcik's shadow server) and SolverMAM.
%       Floor: hol-gap.
%   (c) PRS-relaxing: FCFS preemptive-resume. SolverMVA's PRIOMVA arm and
%       SolverMAM's MMAPPH1PRPR arm -- the latter is phase-type by construction
%       and so is the one competitor that is NOT blind to service moments
%       beyond the mean. Floor: prs-gap.
%   (d) priority-blind but moment-sensitive: Marie's iterative
%       aggregation-decomposition (SolverMVA method 'marie', pfqn_marie) on a
%       plain FCFS station. ansim's closest intellectual relative -- it too
%       solves a small per-station chain against a frozen environment -- but it
%       cannot see the priority order at all. Floor: fcfs-gap.
% Simulation (LDES, SSA) solves PSPRIO natively and needs no relaxation.

here = fileparts(mfilename('fullpath'));
addpath(here);
run('/data/gcasale/line-dev.git/matlab/lineStart.m');

base = struct();
base.D    = [2, 3];
base.Z    = [4, 6];
base.prio = [0, 1];      % lower value = higher priority
base.Wdps = 100;
base.maxIter = 50;
base.tol  = 1e-6;

laws  = {'exp','exp','exp','hyperexp','hyperexp','hyperexp','hypoexp','hypoexp','hypoexp'};
Nvals = [2 4 8 2 4 8 2 4 8];
scv   = containers.Map({'exp','hyperexp','hypoexp'}, {1.0, 4.0, 0.5});

ldesSamples = [1e4, 1e5];
ssaSamples  = 1e5;
seed        = 23000;

R = {};
for c = 1:numel(laws)
    cfg = base;
    cfg.law = laws{c};
    cfg.N   = [Nvals(c), Nvals(c)];

    fprintf('\n=== %s  SCV=%.1f  N=%d ===\n', cfg.law, scv(cfg.law), Nvals(c));

    row = struct('law', cfg.law, 'scv', scv(cfg.law), 'N', Nvals(c));

    % ---- reference -----------------------------------------------------
    [Qref, tref, okref, msgref] = solve_variant(cfg, 'psprio', @(m) SolverCTMC(m,'verbose',false,'force',true));
    if ~okref
        fprintf('  REFERENCE FAILED: %s\n', msgref);
        row.ref_ok = false; R{end+1} = row; continue;
    end
    row.ref_ok = true; row.Qref = Qref; row.t_ref = tref;
    fprintf('  reference PSPRIO CTMC  [%s]  (%.2fs)\n', fmtq(Qref), tref);

    % ---- discipline gaps: exact CTMC of each relaxation ------------------
    row = record(row, 'dpsgap', cfg, Qref, 'dps', @(m) SolverCTMC(m,'verbose',false,'force',true));
    row = record(row, 'holgap', cfg, Qref, 'hol', @(m) SolverCTMC(m,'verbose',false,'force',true));
    row = record(row, 'prsgap', cfg, Qref, 'prs', @(m) SolverCTMC(m,'verbose',false,'force',true));
    row = record(row, 'fcfsgap', cfg, Qref, 'fcfs', @(m) SolverCTMC(m,'verbose',false,'force',true));

    % ---- (a) DPS-emulating analytical methods ---------------------------
    row = record(row, 'FLDmn',  cfg, Qref, 'dps', @(m) SolverFLD(m,'method','minnormal','verbose',false,'force',true));
    row = record(row, 'MVAdps', cfg, Qref, 'dps', @(m) SolverMVA(m,'verbose',false,'force',true));
    row = record(row, 'NCdps',  cfg, Qref, 'dps', @(m) SolverNC(m,'verbose',false,'force',true));

    % ---- (b) HOL-relaxing analytical methods ----------------------------
    row = record(row, 'MVAcl',     cfg, Qref, 'hol', @(m) SolverMVA(m,'config.np_priority','cl','verbose',false,'force',true));
    row = record(row, 'MVAshadow', cfg, Qref, 'hol', @(m) SolverMVA(m,'config.np_priority','shadow','verbose',false,'force',true));
    row = record(row, 'MAMhol',    cfg, Qref, 'hol', @(m) SolverMAM(m,'verbose',false,'force',true));

    % ---- (c) PRS-relaxing analytical methods ----------------------------
    row = record(row, 'MVAprs', cfg, Qref, 'prs', @(m) SolverMVA(m,'verbose',false,'force',true));
    row = record(row, 'MAMprs', cfg, Qref, 'prs', @(m) SolverMAM(m,'verbose',false,'force',true));

    % ---- (d) priority-blind, moment-sensitive ---------------------------
    row = record(row, 'MVAmarie', cfg, Qref, 'fcfs', @(m) SolverMVA(m,'method','marie','verbose',false,'force',true));

    % ---- simulation on the true discipline ------------------------------
    for k = 1:numel(ldesSamples)
        nm = sprintf('LDES%g', ldesSamples(k));
        ns = ldesSamples(k);
        row = record(row, nm, cfg, Qref, 'psprio', @(m) SolverLDES(m,'seed',seed,'samples',ns,'verbose',false,'force',true));
    end
    row = record(row, 'SSA', cfg, Qref, 'psprio', @(m) SolverSSA(m,'seed',seed,'samples',ssaSamples,'verbose',false,'force',true));

    % ---- ansim ----------------------------------------------------------
    for mode = {'ctmcP','live','mn'}
        md = mode{1};
        [Q, t, it, ok, msg] = ansim_solve(cfg, md);
        nm = ['ansim_' md];
        if ok
            e = l1(Q, Qref);
            row.(['err_' nm]) = e; row.(['t_' nm]) = t; row.(['Q_' nm]) = Q; row.(['it_' nm]) = it;
            fprintf('  %-12s L1=%8.4f  (%.2fs, %d iters)  [%s]\n', nm, e, t, it, fmtq(Q));
        else
            row.(['err_' nm]) = NaN; row.(['t_' nm]) = t; row.(['msg_' nm]) = msg;
            fprintf('  %-12s FAILED: %s\n', nm, trunc(msg));
        end
    end

    R{end+1} = row;
    % Checkpoint after every configuration: one slow CTMC late in the sweep must
    % not cost the rows that already finished.
    save(fullfile(here,'poc_ansim_prio.mat'), 'R', 'base', 'laws', 'Nvals');
end

out = fullfile(here, 'poc_ansim_prio.mat');
save(out, 'R', 'base', 'laws', 'Nvals');
fprintf('\nsaved %s\n', out);
ansim_report(R, fullfile(here,'poc_ansim_prio.csv'));
end

% ========================================================================

function row = record(row, name, cfg, Qref, variant, solverfun)
[Q, t, ok, msg] = solve_variant(cfg, variant, solverfun);
if ok
    e = l1(Q, Qref);
    row.(['err_' name]) = e; row.(['t_' name]) = t; row.(['Q_' name]) = Q;
    fprintf('  %-12s L1=%8.4f  (%.2fs)  [%s]\n', name, e, t, fmtq(Q));
else
    row.(['err_' name]) = NaN; row.(['t_' name]) = t; row.(['msg_' name]) = msg;
    fprintf('  %-12s FAILED: %s\n', name, trunc(msg));
end
end

function [Q, t, ok, msg] = solve_variant(cfg, variant, solverfun)
Q = []; ok = true; msg = '';
t0 = tic;
try
    m = ansim_models(cfg, variant);
    Q = solverfun(m).getAvgQLen();
    if any(~isfinite(Q(:)))
        ok = false; msg = 'non-finite queue lengths';
    end
catch ME
    ok = false; msg = ME.message;
end
t = toc(t0);
end

function e = l1(Q, Qref)
if isempty(Q) || ~isequal(size(Q), size(Qref))
    e = NaN;
else
    e = sum(abs(Q(:) - Qref(:)));
end
end

function s = fmtq(Q)
s = strjoin(arrayfun(@(x) sprintf('%.3f',x), Q(:).', 'UniformOutput', false), ' ');
end

function s = trunc(s)
s = strrep(s, newline, ' ');
if numel(s) > 90, s = [s(1:90) '...']; end
end
