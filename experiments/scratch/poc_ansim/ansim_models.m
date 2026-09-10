function [model, cls] = ansim_models(cfg, variant)
% ANSIM_MODELS  Build the Delay+Queue two-class priority model in one of the
% discipline variants used by the ansim priority comparison.
%
%   VARIANT is one of
%     'psprio'  -- the model under study: PS inside the most urgent non-empty
%                  priority group (the reference discipline)
%     'dps'     -- priority emulated by DPS weights w_r = W^(maxprio-prio_r);
%                  this is the only shape SolverFLD and SolverNC can take
%     'hol'     -- non-preemptive head-of-line relaxation (SolverMVA, SolverMAM)
%     'prs'     -- FCFS preemptive-resume relaxation (SolverMVA PRIOMVA arm,
%                  SolverMAM MMAPPH1PRPR)
%     'fcfs'    -- priority dropped entirely, plain FCFS. The only shape
%                  SolverMVA's 'marie' method (Marie's iterative
%                  aggregation-decomposition, pfqn_marie) can take. It is
%                  moment-sensitive but priority-blind, which is exactly the
%                  control that separates the two error sources.
%
% CFG carries D, Z, prio, N, law and Wdps.

M = 2;
K = numel(cfg.N);

model = Network(['prio_' variant]);
nd{1} = Delay(model, 'Think');
switch variant
    case 'psprio', nd{2} = Queue(model, 'Q', SchedStrategy.PSPRIO);
    case 'dps',    nd{2} = Queue(model, 'Q', SchedStrategy.DPS);
    case 'hol',    nd{2} = Queue(model, 'Q', SchedStrategy.HOL);
    case 'prs',    nd{2} = Queue(model, 'Q', SchedStrategy.FCFSPRPRIO);
    case 'fcfs',   nd{2} = Queue(model, 'Q', SchedStrategy.FCFS);
    otherwise, error('unknown variant %s', variant);
end

% DPS carries the priority in its weights, so every class is left at prio 0
% there; the other three variants carry it in the class priority field
% (lower value = higher priority, see solver_amvald_forward.m:39).
if any(strcmp(variant, {'dps','fcfs'}))
    clsprio = zeros(1,K);
else
    clsprio = cfg.prio;
end

cls = cell(1,K);
for r = 1:K
    cls{r} = ClosedClass(model, sprintf('Class%d',r), cfg.N(r), nd{1}, clsprio(r));
end

w = ansim_dpsweights(cfg);
for r = 1:K
    nd{1}.setService(cls{r}, Exp(1/cfg.Z(r)));
    if strcmp(variant, 'dps')
        nd{2}.setService(cls{r}, ansim_svc(cfg.D(r), cfg.law), w(r));
    else
        nd{2}.setService(cls{r}, ansim_svc(cfg.D(r), cfg.law));
    end
end

P = model.initRoutingMatrix;
for r = 1:K
    P{r} = Network.serialRouting(nd);
end
model.link(P);
end
