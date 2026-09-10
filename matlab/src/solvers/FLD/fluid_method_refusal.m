function [bool, reason] = fluid_method_refusal(sn, method, options, model)
% [BOOL, REASON] = FLUID_METHOD_REFUSAL(SN, METHOD, OPTIONS, MODEL)
%
% @brief Every structural rule of a fluid method that the feature registry
% cannot name, asked once for the run and once for the report.
%
% A feature set can refuse a model only for HAVING a registered feature. The
% rules below are of the other kinds -- a station count, a class count, a
% server count, "requires a cache", "requires a patience law", a finite
% horizon, a fork-join route, a binding buffer -- and each of them used to live
% only inside an analyzer arm, where a report could not see it: SolverAUTO's
% findSolver called 'fluid.mol' runnable on a Source -> Delay -> Sink model and
% 'fluid.dae' runnable on an Erlang-fed queue, and each run then stopped.
%
% This function delegates to the per-family predicates, so that each rule has
% one home and two callers: @SolverFLD/runAnalyzer stops on it, and
% @SolverFLD/supportsModelMethod reports it.
%
%   single-station limits   FLUID_QSYS_ADMITS, then FLUID_QSYS_HORIZON for the
%                           three that report a trajectory
%   'mfq' (butools, aoi)    no refusal: off the shape FLUID_MFQ_ADMITS decides
%                           the name RESOLVES to 'matrix' (SolverFLD.resolveMethod),
%                           the documented fallback the analyzer arm keeps; the
%                           predicate is read here only for the AoI exemption
%   'rmf'                   no refusal either: with no Cache node it resolves
%                           to 'matrix', the decomposition's own network step
%   'minnormal', 'refined'  FLUID_MINNORMAL_APPLICABLE (the guards of
%                           SOLVER_FLUID_MOMENTS and FLUID_MOMENT_TERMS)
%   'dae'                   FLUID_PETRI_APPLICABLE on a Petri net, else
%                           FLUID_DAE_APPLICABLE
%   'diffusion'             single or infinite servers only, as
%                           SOLVER_FLUID_DIFFUSION requires
%   every method            FLUID_FORKJOIN_ADMITS, then the binding-capacity
%                           gate (NetworkSolver.checkBindingCapacity)
%
% THE CAPACITY GATE HAS THREE EXEMPTIONS, each because the buffer IS the model
% rather than a limit the drift ignores: 'dae' carries it as an algebraic
% constraint on the drift; 'mol' is stated for the Mt/G/s/0 loss system, where
% the server count is the buffer; and the age-of-information arm of 'mfq' is a
% bufferless or single-buffer queue by definition (AOI_IS_AOI), which the gate
% used to refuse before the arm could run. Nothing else in the FLD tree reads
% sn.cap or sn.classcap.
%
% @param sn NetworkStruct, after sn_nonmarkov_toph, as the run sees it
% @param method the method name, any spelling SolverFLD.canonicalMethod takes
% @param options solver options (horizon, closure limits, rate schedules)
% @param model the Network, optional; only the capacity gate reads it, because
%        it tests the caps the USER set rather than the derived sn.classcap
% @return bool true when METHOD may run on this model
% @return reason the refusal, or '' when BOOL is true

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

bool = true;
reason = '';
m = SolverFLD.canonicalMethod(method);
isAoI = false;
switch m
    case {'ggisgi.fluid','ggingi.tga','tvms','mtginf','mol'}
        [bool, reason] = fluid_qsys_admits(sn, m);
        if bool && any(strcmp(m, {'tvms','mtginf','mol'}))
            [bool, reason] = fluid_qsys_horizon(options);
            if ~bool
                reason = sprintf('The ''%s'' method reports a trajectory. %s', m, reason);
            end
        end
    case 'mfq'
        % A RESOLUTION, NOT A REFUSAL. Off the shape FLUID_MFQ_ADMITS decides,
        % SolverFLD.resolveMethod maps 'mfq' onto 'matrix' (the documented
        % fallback the analyzer arm keeps), so this arm is reached only on a
        % shape it admits, and reads whether it is the age-of-information one.
        % 'rmf' resolves the same way with no Cache node and has no arm here.
        [~, ~, isAoI] = fluid_mfq_admits(sn);
    case {'minnormal','refined'}
        [bool, why] = fluid_minnormal_applicable(sn, options);
        if ~bool
            reason = sprintf('The ''%s'' moment closure cannot answer this model: %s.', m, why);
        end
    case 'dae'
        if any(sn.nodetype == NodeType.Transition)
            [bool, why] = fluid_petri_applicable(sn, options);
        else
            [bool, why] = fluid_dae_applicable(sn, options);
        end
        if ~bool
            reason = sprintf('The ''dae'' method cannot answer this model: %s.', why);
        end
    case 'diffusion'
        multi = find(sn.nservers > 1 & isfinite(sn.nservers), 1);
        if ~isempty(multi)
            bool = false;
            reason = sprintf(['The ''diffusion'' method only supports single-server (c=1) or ' ...
                'infinite-server stations, and station %d has %d servers.'], multi, sn.nservers(multi));
        end
end
if ~bool
    return
end
[bool, reason] = fluid_forkjoin_admits(sn, m);
if ~bool
    return
end
if nargin >= 4 && isa(model, 'Network') && ~any(strcmp(m, {'dae','mol'})) && ~isAoI
    [bool, reason] = NetworkSolver.checkBindingCapacity(model, 'SolverFLD');
    if ~bool
        reason = sprintf(['%s Use options.method=''dae'', which carries the buffer as an ' ...
            'algebraic constraint on the drift.'], reason);
    end
end
end
