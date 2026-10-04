function tf = mva_carries_interlock(sn, options)
% MVA_CARRIES_INTERLOCK True when the MVA path this model already dispatches to
% carries a class-level interlock matrix (Franks 1999, Eq. 4.7) itself, so that
% supplying one does not silently move the model to a DIFFERENT algorithm.
%
% Only two kernels implement the correction: PFQN_MVA (exact, closed
% single-server) and the AMVA forward step of SOLVER_AMVALD. A model that would
% otherwise be solved by exact multiserver or mixed MVA, or by the product-form
% AMVA kernels (linearizer and relatives), cannot take the matrix without
% swapping its algorithm, and the swap is worth far more than the correction it
% carries: inside SolverLN it can turn a converging Picard iteration into a
% limit cycle. A caller holding a matrix such a model cannot carry must apply
% its own correction instead (SolverLN keeps the residt scaling there).
%
% The test is deliberately conservative: where the AMVA product-form branch is
% entered only for some resolved methods, this reports it entered for all of
% them, since refusing the matrix falls back on the caller's own handling
% rather than on a changed algorithm.
%
% Copyright (c) 2012-2026, QORE Lab, Imperial College London
% All rights reserved.

method = 'default';
if isfield(options,'method') && ~isempty(options.method)
    method = regexprep(options.method, '^amva\.', '');
end

hasOpenClass = any(isinf(sn.njobs));
hasClosedClass = any(isfinite(sn.njobs) & sn.njobs > 0);
finiteServers = sn.nservers(isfinite(sn.nservers));
hasFiniteServer = ~isempty(finiteServers);
maxFiniteServers = max(finiteServers); % [] when every station is infinite-server
closedPops = sn.njobs(isfinite(sn.njobs));
closedPopsIntegral = all(closedPops == floor(closedPops));

% pfqn_mva takes the matrix for a closed single-server model, and for nothing else
pfqnMvaCanTakeIt = ~hasOpenClass && closedPopsIntegral && ...
    (isempty(maxFiniteServers) || maxFiniteServers <= 1);

amvaMethods = {'amva','bs','qd','qli','fli','lin','qdlin','sqni','gflin','egflin', ...
    'ab','schmidt','schmidt-ext','tay','scat','aql','qsa', ...
    'lcp','chow','pamb','pami','pamt','clust','dmlin'};

if any(strcmpi(method, {'exact','mva'}))
    tf = pfqnMvaCanTakeIt;
elseif any(strcmpi(method, amvaMethods))
    % Only solver_amvald carries the correction among the AMVA handlers
    tf = ~amva_uses_pf_kernels(sn);
elseif strcmpi(method, 'default')
    if mva_is_bas_model(sn)
        tf = false; % solver_sqd has no interlock term
        return
    end
    exactMixed = hasOpenClass && hasClosedClass && hasFiniteServer && ...
        maxFiniteServers == 1 && sn_has_product_form(sn) && closedPopsIntegral;
    exactSmall = sn.nchains <= 4 && sum(sn.njobs) <= 20 && ...
        sn_has_product_form(sn) && ~sn_has_fractional_populations(sn);
    if exactMixed || exactSmall
        tf = pfqnMvaCanTakeIt;
    else
        tf = ~amva_uses_pf_kernels(sn);
    end
else
    % mvac, sqd, sum, qna, rqna and rqt reach neither kernel
    tf = false;
end
end

function tf = amva_uses_pf_kernels(sn)
% Conservative form of the branch test in SOLVER_AMVA: true when this model may
% be solved by the product-form AMVA kernels rather than by SOLVER_AMVALD. The
% mixed case is reported true for every resolved method, while the branch itself
% takes it only for 'lin', so a caller is never told that solver_amvald will run
% when it might not.
tf = sn_has_product_form_not_het_fcfs(sn) && ...
    ~sn_has_load_dependence(sn) && ...
    isempty(sn.cdscaling) && isempty(sn.jdscaling) && ...
    (~sn_has_open_classes(sn) || sn_has_product_form(sn));
end
