function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_mvald_analyzer(sn, options)
% [Q,U,R,T,C,X,RUNTIME,ITER] = SOLVER_MVALD_ANALYZER(SN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;
method = options.method;

line_debug('MVA load-dependent analyzer starting: method=%s, nclasses=%d, njobs=%s', method, sn.nclasses, mat2str(sn.njobs));

switch options.method
    case {'exact','mva'}
        if ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
            line_error(mfilename,'Exact class-/joint-dependent solver not available in MVA.');
        end
        line_debug('Using exact load-dependent MVA, calling solver_mvald');
        [Q,U,R,T,C,X,lG,iter] = solver_mvald(sn, options);
    case {'default','amva','qd', 'lin', 'qdlin'} %,'aql','qdaql'
        if strcmpi(options.method, 'default') && isempty(sn.cdscaling) && isempty(sn.jdscaling) ...
                && all(isfinite(sn.njobs)) && sn.nchains <= 4 && sum(sn.njobs) <= 20 ...
                && sn_has_product_form(sn) && ~sn_has_fractional_populations(sn)
            % Small closed product-form LD model: use the exact load-dependent
            % MVA recursion, mirroring the default-to-exact upgrade applied to
            % non-LD models in solver_mva_analyzer.
            line_debug('Default method: using exact load-dependent MVA\n');
            [Q,U,R,T,C,X,lG,iter] = solver_mvald(sn, options);
            method = 'exact';
        else
            if strcmpi(options.method, 'default')
                line_debug('Default method: using approximate MVA\n');
            end
            line_debug('Using approximate MVA method: %s, calling solver_amva', options.method);
            [Q,U,R,T,C,X,lG,iter,method] = solver_amva(sn, options);
        end
    otherwise
        line_error(mfilename,sprintf('The %s method is not supported by the load-dependent MVA solver.',options.method));
end

runtime = toc(Tstart);

if options.verbose
    %    line_printf('\nMVA load-dependent analysis completed. Runtime: %f seconds.\n',runtime);
end
return
end
