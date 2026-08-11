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
        elseif isBasModel(sn)
            line_debug('Default method: using BAS for blocking-after-service model\n');
            line_debug(options, 'Model has Blocking-After-Service finite buffers, calling solver_sqd');
            [Q,U,R,T,C,X,lG,iter] = solver_sqd(sn, options);
            method = 'sqd';
        % Force AMVA for class- or joint-dependent models - exact MVA doesn't support them
        elseif ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
            line_debug('Default method: using AMVA for class-/joint-dependent model\n');
            line_debug(options, 'Model has class/joint dependence, calling solver_amva');
            [Q,U,R,T,C,X,lG,iter,method,converged] = solver_amva(sn, options);
        elseif any(isinf(sn.njobs)) && any(isfinite(sn.njobs) & sn.njobs > 0) && max(sn.nservers(isfinite(sn.nservers))) == 1 && sn_has_product_form(sn) && all(sn.njobs(isfinite(sn.njobs)) == floor(sn.njobs(isfinite(sn.njobs))))
            % see _kb/06-solver-catalog.md (MVA section) for the default-dispatch rules
            line_debug('Default method: using exact mixed MVA\n');
            [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
            method = 'exact';
        elseif sn.nchains <= 4 && sum(sn.njobs) <= 20 && sn_has_product_form(sn) && ~sn_has_fractional_populations(sn)
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
    case {'amva','bs','qd','qli','fli','lin','qdlin','sqni','egflin','gflin','ab','schmidt','schmidt-ext'} %,'aql','qdaql'
        line_debug(options, 'Using approximate MVA method: %s, calling solver_amva', method);
        [Q,U,R,T,C,X,lG,iter,~,converged] = solver_amva(sn, options);
    otherwise
        % An unsupported method must error by name (returning empty metrics let the
        % caller index into [] with an opaque message). 'aql'/'qdaql' exist in
        % solver_amva but are deliberately not dispatched here.
        line_error(mfilename, sprintf(['The ''%s'' method is not dispatched by ' ...
            'solver_mva_analyzer. Supported: default, exact, mva, mvac, amva, bs, qd, qli, fli, ' ...
            'lin, qdlin, sqni, egflin, gflin, ab, schmidt, schmidt-ext.'], method));
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

function tf = isBasModel(sn)
% Detect a closed single-chain network with Blocking-After-Service (BAS) finite-buffer
% blocking, which solver_sqd handles but exact/AMVA MVA does not.
tf = false;
if sn.nchains ~= 1 || sn.nclosedjobs <= 0 || isempty(sn.droprule)
    return;
end
if any(isinf(sn.njobs))
    return; % open class present
end
tf = any(sn.droprule(:) == DropStrategy.BAS);
end
