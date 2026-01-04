function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_mva_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_MVA_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter = NaN;
Tstart = tic;
method = options.method;
method = regexprep(method, '^amva\.', '');

line_debug(options, 'MVA analyzer starting: method=%s, nclasses=%d, njobs=%s, nchains=%d', method, sn.nclasses, mat2str(sn.njobs), sn.nchains);

switch method
    case {'exact','mva'}
        line_debug(options, 'Using exact MVA method, calling solver_mva');
        [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
    case {'qna'}
        line_debug(options, 'Using QNA method, calling solver_qna');
        [Q,U,R,T,C,X] = solver_qna(sn, options);
        lG = NaN;
    case {'default'}
        % for non-exponential open queueing networks, use qna
        % (commented as not ready yet, it fails on example_cacheModel_3.m)
        %if all(isinf(sn.njobs)) && any(any(sn.scv ~= 1.0))
        %    line_warning(mfilename,'QNA implementation is still in beta version.')
        %    [Q,U,R,T,C,X] = solver_qna(sn, options);
        %    lG = NaN;
        %else
        % Force AMVA for joint-dependent (LJD/LJCD) models - exact MVA doesn't support it
        if sn_has_joint_dependence(sn)
            line_debug('Default method: using AMVA for joint-dependent model\n');
            line_debug(options, 'Model has joint dependence (LJD/LJCD), calling solver_amva');
            [Q,U,R,T,C,X,lG,iter,method] = solver_amva(sn, options);
        elseif sn.nchains <= 4 && sum(sn.njobs) <= 20 && sn_has_product_form(sn) && ~sn_has_fractional_populations(sn)
            line_debug('Default method: using exact MVA\n');
            % The parameters above take in the worst case a handful of ms
            line_debug(options, 'Model qualifies for exact MVA (nchains=%d, njobs=%d, product-form=%d), calling solver_mva', sn.nchains, sum(sn.njobs), sn_has_product_form(sn));
            [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
            method = 'exact';
        else
            line_debug('Default method: using approximate MVA\n');
            line_debug(options, 'Model requires approximation (nchains=%d, njobs=%d, product-form=%d), calling solver_amva', sn.nchains, sum(sn.njobs), sn_has_product_form(sn));
            [Q,U,R,T,C,X,lG,iter,method] = solver_amva(sn, options);
        end
        %end
    case {'amva','bs','qd','qli','fli','lin','qdlin','sqni','egflin','gflin','ab','schmidt','schmidt-ext'} %,'aql','qdaql'
        line_debug(options, 'Using approximate MVA method: %s, calling solver_amva', method);
        [Q,U,R,T,C,X,lG,iter] = solver_amva(sn, options);
    otherwise
        Q=[];
        U=[];
        R=[];
        T=[];
        C=[];
        X=[];
        lG=[];
        iter=[];
        if options.verbose
            line_warning(mfilename,'Unsupported SolverMVA method.');
        end
end
runtime = toc(Tstart);

if options.verbose
    %line_printf('\nMVA analysis completed. Runtime: %f seconds.\n',runtime);
end

end
