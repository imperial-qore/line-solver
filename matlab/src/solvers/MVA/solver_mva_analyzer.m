function [Q,U,R,T,C,X,lG,runtime,iter] = solver_mva_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME] = SOLVER_MVA_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

iter = NaN;
Tstart = tic;
switch options.method
    case {'exact','mva'}
        [Q,U,R,T,C,X,lG] = solver_mva(sn, options);
    case {'qna'}
        line_warning(mfilename,'QNA implementation is still in beta version.\n')
        [Q,U,R,T,C,X] = solver_qna(sn, options);
        lG = NaN;
    case {'default'}
        % for non-exponential open queueing networks, use qna
        % (commented as not ready yet, it fails on example_cacheModel_5.m)
        %if all(isinf(sn.njobs)) && any(any(sn.scv ~= 1.0))
        %    line_warning(mfilename,'QNA implementation is still in beta version.')
        %    [Q,U,R,T,C,X] = solver_qna(sn, options);
        %    lG = NaN;
        %else
        [Q,U,R,T,C,X,lG,iter] = solver_amva(sn, options);
        %end
    case {'amva','bs','qd','qli','fli','aql','qdaql','lin','qdlin','sqni','egflin','gflin'}
        [Q,U,R,T,C,X,lG,iter] = solver_amva(sn, options);
    otherwise
        line_error(mfilename,'Unsupported SolverMVA method.');
end
runtime = toc(Tstart);

if options.verbose
    %line_printf('\nMVA analysis completed. Runtime: %f seconds.\n',runtime);
end

end
