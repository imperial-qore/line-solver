function [Q,U,R,T,C,X,lG,runtime,totiter,actualmethod] = solver_mva_sjn_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,TOTITER,ACTUALMETHOD] = SOLVER_MVA_SJN_ANALYZER(SN, OPTIONS)
%
% Analyzer for closed networks with non-preemptive shortest-job-next (SJF)
% stations. The station is modelled by the conditional waiting time equation
% of K. Kant, "MVA approximations for SJN scheduling", Performance Evaluation
% 15(1):41-61, 1992, evaluated either over the full population lattice
% (pfqn_mvasjn) or through its Schweitzer fixed point (pfqn_amvasjn).
%
% The lattice costs prod(N+1) steps, so 'default' switches to the fixed point
% once the lattice exceeds options.config.sjn_lattice_max states. Ask for
% 'exact' or 'mva' to force the lattice, 'amva' to force the fixed point.
% 'default' also falls back to the fixed point when the lattice reports the
% starvation regime, where it has no valid continuation.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
[Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = sn_get_demands_chain(sn); %#ok<ASGLU>

M = sn.nstations;
K = sn.nchains;
sched = sn.sched;
nservers = sn.nservers;

% The premise (closed model, single-server stations with BCMP disciplines, the
% method names this analyzer dispatches) is judged by SolverMVA.supportsSjn, the
% predicate the report gates on and mvaDispatch refuses by, so the sentence a
% caller reads here is the sentence that kept the pair off the report.
[sjnOk, sjnReason] = SolverMVA.supportsSjn(sn, options.method);
if ~sjnOk
    line_error(mfilename, sjnReason);
end

rows = [];      % queueing stations passed to the SJN recursion
infrows = [];   % delay stations, folded into the think time
sjnrows = [];   % positions within rows of the SJN stations
for ist = 1:M
    switch sched(ist)
        case SchedStrategy.EXT
            % no-op, there is no external arrival stream in a closed model
        case SchedStrategy.INF
            infrows(1,end+1) = ist; %#ok<AGROW>
        case SchedStrategy.SJF
            rows(1,end+1) = ist; %#ok<AGROW>
            sjnrows(1,end+1) = length(rows); %#ok<AGROW>
        otherwise % PS, LCFSPR, FCFS, SIRO: anything else was refused above
            rows(1,end+1) = ist; %#ok<AGROW>
    end
end

L = STchain(rows,:) .* Vchain(rows,:);
V = Vchain(rows,:);
scv = ones(length(rows),K);
for j = sjnrows
    ist = rows(j);
    for r = 1:K
        if isfinite(SCVchain(ist,r)) && SCVchain(ist,r) > 0
            scv(j,r) = SCVchain(ist,r);
        end
    end
end
Z = zeros(1,K);
for j = 1:length(infrows)
    Z = Z + STchain(infrows(j),:) .* Vchain(infrows(j),:);
end

sjnopt = struct();
if isfield(options,'iter_tol') && ~isempty(options.iter_tol)
    sjnopt.tol = options.iter_tol;
end
if isfield(options,'iter_max') && ~isempty(options.iter_max)
    sjnopt.iter_max = options.iter_max;
end
if isfield(options,'config')
    if isfield(options.config,'sjn_ns'), sjnopt.ns = options.config.sjn_ns; end
    if isfield(options.config,'sjn_lfactor'), sjnopt.Lfactor = options.config.sjn_lfactor; end
    if isfield(options.config,'sjn_umax'), sjnopt.umax = options.config.sjn_umax; end
end
% SJN applies within a class and the classes are then non-preemptively prioritised;
% without distinct priorities the jobs of every class are compared by size directly
if sn.nchains == sn.nclasses && length(unique(sn.classprio)) == sn.nclasses
    sjnopt.prio = sn.classprio(:)';
end

latticemax = 1e5;
if isfield(options,'config') && isfield(options.config,'sjn_lattice_max')
    latticemax = options.config.sjn_lattice_max;
end
switch options.method
    case {'amva','bs','sjn.amva'}
        uselattice = false;
    case {'exact','mva','sjn.mva'}
        uselattice = true;
    otherwise
        uselattice = prod(Nchain+1) <= latticemax;
end

if uselattice
    try
        [Xchain,Qrows,Urows] = pfqn_mvasjn(L,Nchain,Z,scv,sjnrows,V,sjnopt);
        totiter = 1;
        actualmethod = 'sjn.mva';
    catch ME
        if ~strcmp(ME.identifier,'LINE:SjnStarvation') || ~strcmp(options.method,'default')
            rethrow(ME);
        end
        line_debug(options, 'SJN station in the starvation regime, falling back to the fixed point');
        [Xchain,Qrows,Urows,~,~,totiter] = pfqn_amvasjn(L,Nchain,Z,scv,sjnrows,V,sjnopt);
        actualmethod = 'sjn.amva';
    end
else
    [Xchain,Qrows,Urows,~,~,totiter] = pfqn_amvasjn(L,Nchain,Z,scv,sjnrows,V,sjnopt);
    actualmethod = 'sjn.amva';
end

Qchain = zeros(M,K);
Uchain = zeros(M,K);
Tchain = zeros(M,K);
Qchain(rows,:) = Qrows;
Uchain(rows,:) = Urows;
for j = 1:length(infrows)
    ist = infrows(j);
    Qchain(ist,:) = Xchain .* STchain(ist,:) .* Vchain(ist,:);
    Uchain(ist,:) = Qchain(ist,:);
end
for k = 1:M
    for r = 1:K
        Tchain(k,r) = Xchain(r) * Vchain(k,r);
    end
end
Rchain = Qchain ./ Tchain;

Xchain(~isfinite(Xchain)) = 0;
Uchain(~isfinite(Uchain)) = 0;
Qchain(~isfinite(Qchain)) = 0;
Rchain(~isfinite(Rchain)) = 0;
Xchain(Nchain==0) = 0;
Uchain(:,Nchain==0) = 0;
Qchain(:,Nchain==0) = 0;
Rchain(:,Nchain==0) = 0;
Tchain(:,Nchain==0) = 0;

lG = NaN;
[Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
runtime = toc(T0);
end
