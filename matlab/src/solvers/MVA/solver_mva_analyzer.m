function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_mva_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_MVA_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter = NaN;
% Convergence flag as reported by the handler; [] means it reports none,
% in which case the count-vs-budget test below applies instead.
converged = [];
Tstart = tic;
method = options.method;
method = regexprep(method, '^amva\.', '');

line_debug(options, 'MVA analyzer starting: method=%s, nclasses=%d, njobs=%s, nchains=%d', method, sn.nclasses, mat2str(sn.njobs), sn.nchains);

switch method
    case {'exact','mva'}
        line_debug(options, 'Using exact MVA method, calling solver_mva');
        [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
    case {'mvac'}
        line_debug(options, 'Using exact MVAC method, calling solver_mvac');
        [Q,U,R,T,C,X,lG] = solver_mvac(sn, options);
    case {'sqd'}
        line_debug(options, 'Using BAS method, calling solver_sqd');
        [Q,U,R,T,C,X,lG,iter] = solver_sqd(sn, options);
    case {'qna'}
        line_debug(options, 'Using QNA method, calling solver_qna');
        [Q,U,R,T,C,X] = solver_qna(sn, options);
        lG = NaN;
    case {'sum','esum'}
        line_debug(options, 'Using summation method, calling solver_mva_sum');
        [Q,U,R,T,C,X,lG,iter] = solver_mva_sum(sn, options);
    case {'rqna'}
        line_debug(options, 'Using RQNA method, calling solver_rqna');
        [Q,U,R,T,C,X] = solver_rqna(sn, options);
        lG = NaN;
    case {'rqt'}
        line_debug(options, 'Using RQT method, calling solver_rqt');
        [Q,U,R,T,C,X] = solver_rqt(sn, options);
        lG = NaN;
    case {'default'}
        % for non-exponential open queueing networks, use qna
        % (commented as not ready yet, it fails on example_cacheModel_3.m)
        %if all(isinf(sn.njobs)) && any(any(sn.scv ~= 1.0))
        %    line_warning(mfilename,'QNA implementation is still in beta version.')
        %    [Q,U,R,T,C,X] = solver_qna(sn, options);
        %    lG = NaN;
        %else
        % see _kb/06-solver-catalog.md (MVA section) for the default-dispatch rules
        if sn.nclasses == 1 && all(isinf(sn.njobs)) && sn_has_bursty_arrival(sn)
            line_debug('Default method: bursty single-class open network, using RQNA\n');
            line_debug(options, 'Non-renewal arrivals detected, calling solver_rqna');
            [Q,U,R,T,C,X] = solver_rqna(sn, options);
            lG = NaN;
            method = 'rqna';
        % Closed single-chain Blocking-After-Service network (finite buffers)
        elseif mva_is_bas_model(sn)
            line_debug('Default method: using BAS for blocking-after-service model\n');
            line_debug(options, 'Model has Blocking-After-Service finite buffers, calling solver_sqd');
            [Q,U,R,T,C,X,lG,iter] = solver_sqd(sn, options);
            method = 'sqd';
        % An LCFS station has one arm, the exact LCFS/LCFS-PR recursion in
        % solver_mva (SolverMVA.supportsLcfs has admitted the pair); the
        % population thresholds below must not send it to AMVA, which has no
        % LCFS arm and returned a zero wait at the station.
        elseif any(sn.sched == SchedStrategy.LCFS)
            line_debug('Default method: using the exact LCFS/LCFS-PR recursion\n');
            [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
            method = 'exact';
        % Force AMVA for class- or joint-dependent models - exact MVA doesn't support them
        elseif ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
            line_debug('Default method: using AMVA for class-/joint-dependent model\n');
            line_debug(options, 'Model has class/joint dependence, calling solver_amva');
            [Q,U,R,T,C,X,lG,iter,method,converged] = solver_amva(sn, options);
        elseif ~needsAmvaForInterlock(sn, options) && any(isinf(sn.njobs)) && any(isfinite(sn.njobs) & sn.njobs > 0) && max(sn.nservers(isfinite(sn.nservers))) == 1 && sn_has_product_form(sn) && all(sn.njobs(isfinite(sn.njobs)) == floor(sn.njobs(isfinite(sn.njobs))))
            % see _kb/06-solver-catalog.md (MVA section) for the default-dispatch rules
            line_debug('Default method: using exact mixed MVA\n');
            [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
            method = 'exact';
        elseif ~needsAmvaForInterlock(sn, options) && sn.nchains <= 4 && sum(sn.njobs) <= 20 && sn_has_product_form(sn) && ~sn_has_fractional_populations(sn)
            line_debug('Default method: using exact MVA\n');
            % The parameters above take in the worst case a handful of ms
            line_debug(options, 'Model qualifies for exact MVA (nchains=%d, njobs=%d, product-form=%d), calling solver_mva', sn.nchains, sum(sn.njobs), sn_has_product_form(sn));
            [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
            method = 'exact';
        else
            line_debug('Default method: using approximate MVA\n');
            line_debug(options, 'Model requires approximation (nchains=%d, njobs=%d, product-form=%d), calling solver_amva', sn.nchains, sum(sn.njobs), sn_has_product_form(sn));
            [Q,U,R,T,C,X,lG,iter,method,converged] = solver_amva(sn, options);
        end
        %end
    case {'amva','bs','qd','qli','fli','lin','qdlin','sqni','egflin','gflin','ab','schmidt','schmidt-ext','tay','scat','aql','qsa', ...
            'lcp','chow','pamb','pami','pamt','clust','dmlin','priomva'}
        line_debug(options, 'Using approximate MVA method: %s, calling solver_amva', method);
        [Q,U,R,T,C,X,lG,iter,~,converged] = solver_amva(sn, options);
    otherwise
        % An unsupported method must error by name (returning empty metrics let the
        % caller index into [] with an opaque message). 'qdaql' is NOT dispatched:
        % solver_amvald_forward has no aql arm, so the name would silently compute
        % qd (the trap documented at solver_amva.m:54 for amva.tay).
        line_error(mfilename, sprintf(['The ''%s'' method is not dispatched by ' ...
            'solver_mva_analyzer. Supported: default, exact, mva, mvac, amva, aql, qsa, bs, qd, qli, fli, ' ...
            'lin, qdlin, sqni, egflin, gflin, ab, schmidt, schmidt-ext, tay, scat, ' ...
            'lcp, chow, pamb, pami, pamt, clust, dmlin, priomva.'], method));
end
runtime = toc(Tstart);

% see _kb/06-solver-catalog.md (MVA section) for AMVA convergence flag-vs-count
iterBudget = min(options.iter_max, 10000);
if ~isempty(converged)
    nonConverged = ~converged;
else
    nonConverged = ~isempty(iter) && isscalar(iter) && ~isnan(iter) && isfinite(iter) && iter >= iterBudget;
end
if nonConverged
    line_warning_always(mfilename, ...
        'AMVA method ''%s'' did not meet the convergence tolerance %g after %d iterations; the returned metrics may not be converged. Try another method (e.g. ''qd'' or ''bs''), raise options.iter_max, or loosen options.iter_tol.', ...
        method, options.iter_tol, iter);
end

if options.verbose
    %line_printf('\nMVA analysis completed. Runtime: %f seconds.\n',runtime);
end

end

function tf = needsAmvaForInterlock(sn, options)
% True when an interlock matrix is supplied but exact MVA cannot honour it.
% PFQN_MVA carries the Eq. (4.7) correction for closed single-server models
% only, so a mixed or multiserver model with an interlock goes to AMVA, which
% applies the same correction to the arrival-instant queue length.
tf = false;
if ~isfield(options,'config') || ~isfield(options.config,'interlock') || isempty(options.config.interlock)
    return
end
tf = any(isinf(sn.njobs)) || max(sn.nservers(isfinite(sn.nservers))) > 1;
end
