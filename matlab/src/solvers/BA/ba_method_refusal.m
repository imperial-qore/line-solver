function reason = ba_method_refusal(sn, method, options)
% REASON = BA_METHOD_REFUSAL(SN, METHOD, OPTIONS)
%
% The STRUCTURAL premises of the SolverBA bound families, in one place: the
% reason METHOD cannot bound the model SN, or '' when it can.
%
% ONE PREDICATE, TWO CALLERS. @SolverBA/runAnalyzer and solver_ba_analyzer ask
% it before dispatching and raise the reason it returns; SolverBA.supportsModelMethod
% asks it after the feature gate and reports the same sentence, which is what
% findSolver, listValidMethods and SolverAUTO's ranked choice all read. A second
% copy of any rule below is how the report and the run drift apart: the report
% offers a pair that raises the moment it is run, which is the defect this
% function exists to remove.
%
% WHAT BELONGS HERE AND WHAT DOES NOT. Only the rules the feature registry
% cannot name. SolverFeatureSet.fields has no entry for "one class", for a
% server count, for a station count, for a binding buffer, for "FCFS with a
% non-exponential law" (a conjunction of two names the registry holds
% separately) or for "a Transition node is present", so those rules are
% structural and live here. Rules of the form "this family does not accept a
% delay station", "does not accept an open class" or "does not read a MAP" ARE
% nameable, and belong in SolverBA.getMethodFeatureSet, which drops
% SchedStrategy_INF / ClosedClass / MAP from the offending method's set instead
% -- a feature set can refuse a model for HAVING a construct, never for
% lacking one.
%
% THE PREMISES, family by family:
%   single-class closed  aba, bjb, pb, sb, gb, harel, lr, pbh, cbh, pbk, bjbk,
%                        ssd, sib, scb, ldbcmp, mapamva, qr/qrf.* and the auto
%                        composite. Every one of them is a function of the
%                        single-chain demand vector D = V./rates (or of the
%                        single-chain phase-type chain, for the QRF and
%                        MAP-AMVA reductions), the think time Z and the
%                        population N, which a multiclass or open model does
%                        not have.
%   fully closed         mwba, cub, mbjb, looping. These are the multiclass
%                        families: they take a per-chain demand matrix and a
%                        population VECTOR, so several classes are fine but an
%                        infinite population is not.
%   fully open           bpt, bgt, snc, judged in BA_OPEN_REFUSAL together
%                        with the routing premises of bgt and snc.
%   one server           every family above except ssd, ldbcmp, auto and the
%                        two load-dependent QRF arms. ssd is the multiserver
%                        bound (Suri-Dallery), ldbcmp is a lower bound that
%                        stays one when a station gains servers, auto composes
%                        whichever candidates survive, and 'qrf.mmi.ld' /
%                        'qrf.mmi.linear' carry alpha(i,n) = min(n,c).
%   product form         every family parameterized by demands -- all of the
%                        single-class closed and fully closed ones above except
%                        mwba and the phase-type reductions -- is derived for a
%                        BCMP network and reads the service MEAN alone, so a
%                        FCFS station whose law is not exponential is refused:
%                        the arm would return the bound of the exponential
%                        network, which is a bracket for a different system. PS,
%                        LCFSPR and INF are insensitive and admit any law. The
%                        three multiclass BCMP families (cub, mbjb, looping)
%                        also refuse class-dependent FCFS rates, which BCMP
%                        type 1 does not allow. mwba is the exception: it is
%                        derived for NBUE service at FIFO and priority stations,
%                        so there it admits Exp, Erlang, Det and Uniform and
%                        refuses a law the registry cannot certify as NBUE.
%   phase-type chain     qr/qrf.* and mapamva. Their chain carries ONE phase
%                        per station, the phase of the job in service, which
%                        is FCFS semantics: a multi-phase law at a PS or LCFSPR
%                        single server would describe a different chain (that
%                        station's true law is the exponential one, by
%                        insensitivity), so it is refused. SN_TO_QRF_ALPHA owns
%                        the delay / multiserver / load-dependent rate law of
%                        the QRF arms and the one restriction that survives
%                        there (exponential service where a station serves
%                        several jobs at once); the alpha-free arms refuse any
%                        such station by naming the two arms that serve it. The
%                        blocking arms need the tables SN_TO_QRF_BLOCKING
%                        derives and the capacities SN_TO_QRF_CAPACITY reads.
%                        mapamva additionally carries phases at ONE station.
%   finite buffer        every family but the QRF blocking bounds and spnlp
%                        (BA_IGNORES_BLOCKING) presumes unbounded buffers, so a
%                        buffer that BINDS refuses it. 'default', 'auto' and
%                        'auto.upper' are first routed to 'qrf.bas' where that
%                        bound applies (BA_RESOLVE_MODEL_METHOD), and the
%                        refusal otherwise says why the routing did not.
%   Petri net            spnlp.* alone bounds a model holding Transition nodes,
%                        and bounds nothing else; its own mode rules are in
%                        BA_SPNLP_REFUSAL.
%
% METHOD is taken as the caller spells it and resolved through
% BA_RESOLVE_MODEL_METHOD, so 'default' is judged as the gb.upper it runs as
% (or as the qrf.bas it runs as on a blocked model) and the reason names that.
% OPTIONS is optional and is read for the QRF overrides (config.qrf_params,
% config.qrf_maxvars) exactly as solver_ba_qrf_analyzer reads them.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

reason = '';
if nargin < 3
    options = [];
end
[resolved, blockWhy] = ba_resolve_model_method(sn, method);
dot = strfind(resolved,'.');
if isempty(dot)
    fam = resolved;
else
    fam = resolved(1:dot(1)-1);
end

% -- the Petri-net split, total in both directions ----------------------
if strcmp(fam, 'spnlp')
    reason = ba_spnlp_refusal(sn, resolved);
    return
elseif any(sn.nodetype == NodeType.Transition)
    reason = sprintf(['Method ''%s'' is parameterized by demands and a population, which the ' ...
        'marking of a Petri net is not. Use ''spnlp.upper''/''spnlp.lower'' or the ' ...
        '''spnlp.op.*'' pair on this model.'], resolved);
    return
end

% -- finite-buffer blocking ----------------------------------------------
% Refusing is not conservatism: on cqn_bas_blocking (Queue2 capped at 1,
% N = 2) gb.upper reports QLen 1.28 at a station that can never hold more than
% one job, the unblocked model's answer verbatim.
if ba_ignores_blocking(resolved) && sn_has_blocking(sn)
    if isempty(blockWhy)
        detail = '';
    else
        detail = sprintf(' The QRF blocking bounds do not apply here either: %s', blockWhy);
    end
    reason = sprintf(['Method ''%s'' does not support finite-buffer blocking: every SolverBA ' ...
        'bound family but the QRF blocking ones is parameterized by demands and a population ' ...
        'alone, so it bounds the model as if its buffers were unbounded. Use SolverMVA with ' ...
        'method ''sqd'', an exact solver (CTMC, SSA, JMT, LDES), or the QRF blocking bounds ' ...
        '''qrf.bas''/''qrf.rsrd'', which model the finite buffer.%s'], resolved, detail);
    return
end

singleClassFams = {'auto','aba','bjb','pb','sb','gb','harel','lr', ...
    'pbh','cbh','pbk','bjbk','ssd','sib','scb','ldbcmp','mapamva','qrf'};
fullyClosedFams = {'mwba','cub','mbjb','looping'};
% The families whose alternative on a multiserver model IS 'ssd', which is why
% the reason names it. The multiclass ones below have no such counterpart, and
% the QRF arms name their own load-dependent siblings instead.
singleServerFams = {'aba','bjb','pb','sb','gb','harel','lr', ...
    'pbh','cbh','pbk','bjbk','sib','scb','mapamva'};
% The BCMP-derived families, which read the service mean alone.
productFormFams = {'auto','aba','bjb','pb','sb','gb','harel','lr', ...
    'pbh','cbh','pbk','bjbk','ssd','sib','scb','ldbcmp','cub','mbjb','looping'};

if any(strcmp(fam, singleClassFams))
    if sn.nclasses ~= 1 || sn.nclosedjobs <= 0
        reason = sprintf('Method ''%s'' supports single-class closed networks only.', resolved);
        return
    end
elseif any(strcmp(fam, fullyClosedFams))
    if sn.nclosedjobs <= 0 || any(isinf(sn.njobs))
        reason = sprintf('Method ''%s'' supports fully closed networks only.', resolved);
        return
    end
elseif any(strcmp(fam, {'bpt','bgt','snc'}))
    reason = ba_open_refusal(sn, resolved);
    return
end

if any(sn.nservers(sn.sched ~= SchedStrategy.INF) > 1)
    if any(strcmp(fam, singleServerFams))
        reason = sprintf(['Method ''%s'' does not support multi-server stations ' ...
            '(use ''ssd'').'], resolved);
        return
    elseif any(strcmp(fam, fullyClosedFams))
        reason = sprintf('Method ''%s'' does not support multi-server stations.', resolved);
        return
    end
end

% -- product form: FCFS needs exponential service --------------------------
% Judged over the pairs a FCFS station could serve (finite positive rate), the
% same outer approximation SN_HAS_PRODUCT_FORM applies, and by ProcessType
% rather than by SCV so that the refusal names the law the model declares.
if any(strcmp(fam, productFormFams))
    for i = find(sn.sched(:)' == SchedStrategy.FCFS)
        for r = 1:sn.nclasses
            if ~isfinite(sn.rates(i,r)) || sn.rates(i,r) <= 0
                continue
            end
            if sn.procid(i,r) ~= ProcessType.EXP
                reason = sprintf(['Method ''%s'' is derived for a product-form network and reads ' ...
                    'the service mean alone: station %d class %d serves FCFS with %s service, so ' ...
                    'the bound would bracket the exponential network instead. Use PS or LCFSPR at ' ...
                    'that station, exponential service, or the QRF family (''qr'', ''qrf.mmi''), ' ...
                    'whose chain carries the phase-type law.'], resolved, i, r, ...
                    ProcessType.toText(sn.procid(i,r)));
                return
            end
        end
    end
    if any(strcmp(fam, {'cub','mbjb','looping'})) && sn_has_multi_class_heter_fcfs(sn)
        reason = sprintf(['Method ''%s'' is derived for a product-form multiclass network, ' ...
            'and BCMP type 1 requires a FCFS station to serve every class at one rate; this ' ...
            'model has class-dependent FCFS rates. Use PS at that station, or ''mwba'', which ' ...
            'is derived for class-dependent FIFO demands.'], resolved);
        return
    end
end

% -- mwba: NBUE service wherever the discipline is not insensitive ---------
% Majumdar-Woodside assume NBUE service; PS, LCFSPR and INF are insensitive
% (their exact solution is the exponential one, which is NBUE), everywhere
% else the declared law has to be NBUE by itself.
if strcmp(fam, 'mwba')
    nbue = [ProcessType.EXP, ProcessType.ERLANG, ProcessType.DET, ProcessType.UNIFORM];
    insensitive = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.LCFSPR];
    for i = 1:sn.nstations
        if any(sn.sched(i) == insensitive)
            continue
        end
        for r = 1:sn.nclasses
            if ~isfinite(sn.rates(i,r)) || sn.rates(i,r) <= 0
                continue
            end
            if ~any(sn.procid(i,r) == nbue)
                reason = sprintf(['Method ''%s'' assumes NBUE service at a station that is not ' ...
                    'insensitive: station %d class %d serves %s with %s service, which the ' ...
                    'registry cannot certify as NBUE. Use Exp, Erlang, Det or Uniform there, or ' ...
                    'PS/LCFSPR scheduling.'], resolved, i, r, SchedStrategy.toText(sn.sched(i)), ...
                    ProcessType.toText(sn.procid(i,r)));
                return
            end
        end
    end
end

% -- the phase-type chains: qr/qrf.* and mapamva ---------------------------
if strcmp(fam, 'qrf') || strcmp(fam, 'mapamva')
    if strcmp(fam, 'qrf')
        [~, alphaMsg, isLd] = sn_to_qrf_alpha(sn);
        ldArm = any(strcmp(resolved, {'qrf.mmi.ld','qrf.mmi.linear'}));
        if isLd && ~ldArm
            reason = sprintf(['the ''%s'' method models every station as a single server: its ' ...
                'transition rates carry no population index, so it has nowhere to put the rate ' ...
                'of a delay, a multiserver or a load-dependent station. Use ''qrf.mmi.ld'' or ' ...
                '''qrf.mmi.linear'', which do.'], resolved);
            return
        end
        if ~isempty(alphaMsg)
            reason = sprintf('The ''%s'' method cannot be applied: %s', resolved, alphaMsg);
            return
        end
    end
    % One phase per station is the phase of the job in service, i.e. FCFS. A
    % delay or a multiserver with phases is refused above by sn_to_qrf_alpha for
    % the QRF arms, and mapamva refuses a delay through its feature set.
    phased = [];
    for i = 1:sn.nstations
        ki = 1;
        if ~isempty(sn.proc) && numel(sn.proc) >= i && ~isempty(sn.proc{i}) && ~isempty(sn.proc{i}{1})
            ki = size(sn.proc{i}{1}{1}, 1);
        end
        if ki > 1
            phased(end+1) = i; %#ok<AGROW>
            if sn.sched(i) ~= SchedStrategy.FCFS && sn.sched(i) ~= SchedStrategy.INF && ...
                    sn.nservers(i) <= 1
                reason = sprintf(['Method ''%s'' carries the phase of the job in service at each ' ...
                    'station, which is FCFS semantics: station %d serves %s with a %d-phase law, ' ...
                    'whose exact solution is the exponential one by insensitivity and not the ' ...
                    'chain this bound relaxes. Use FCFS at that station, or exponential service.'], ...
                    resolved, i, SchedStrategy.toText(sn.sched(i)), ki);
                return
            end
        end
    end
    if strcmp(fam, 'mapamva') && numel(phased) > 1
        reason = sprintf(['Method ''%s'' carries phases at ONE station: the LP gives queue M the ' ...
            '(D0,D1) pair and every other queue a scalar rate. Stations %s are all ' ...
            'non-exponential. Use a QRF method, whose q carries a phase at every station.'], ...
            resolved, mat2str(phased(:)'));
        return
    end
    if startsWith(resolved, 'qrf.bas') || strcmp(resolved, 'qrf.rsrd')
        [~, ~, capMsg] = sn_to_qrf_capacity(sn);
        if ~isempty(capMsg)
            reason = sprintf('The ''%s'' method cannot be applied: %s', resolved, capMsg);
            return
        end
    end
    if startsWith(resolved, 'qrf.bas')
        % The tables are derived from the model unless the caller supplies
        % them, the override solver_ba_qrf_analyzer honours as well.
        supplied = isstruct(options) && isfield(options,'config') && isstruct(options.config) && ...
            isfield(options.config,'qrf_params') && ~isempty(options.config.qrf_params);
        if ~supplied
            [~, blkMsg] = sn_to_qrf_blocking(sn, options);
            if ~isempty(blkMsg)
                reason = sprintf(['The ''%s'' method cannot be applied to this model: %s ' ...
                    'Supply options.config.qrf_params explicitly to override the derivation.'], ...
                    resolved, blkMsg);
                return
            end
        end
    end
end
end
