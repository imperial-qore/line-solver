function [tf, reason] = mam_retrial_applicable(sn)
% [TF, REASON] = MAM_RETRIAL_APPLICABLE(SN)
%
% Can the 'retrial' method of SolverMAM answer this model? REASON is '' when it
% can, and otherwise names what is missing.
%
% THE RULE IS A "MUST BE PRESENT" ONE, which is why it cannot live in a feature
% set: a SolverFeatureSet says "I accept this construct", so it can refuse a
% model for HAVING something and never for LACKING it. solver_mam_retrial needs
% an impatience configuration to analyze -- either the BMAP/PH/N/N bufferless
% retrial topology of Dudin et al. (Mathematics 13(9), 2025) or a reneging
% patience law for the MAP/M/s+G analysis of Gursoy, Mehr and Akar -- and a
% model carrying neither is not a smaller retrial model, it is a different one.
%
% ONE PREDICATE, TWO CALLERS: SolverMAM.supportsModelMethod asks it so the
% method is not offered by findSolver or SolverAUTO on a model it cannot
% answer, and solver_mam_analyzer asks it so a caller naming 'retrial' by hand
% gets the identical sentence.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[isRetrial, retInfo] = qsys_is_retrial(sn);
% The reneging branch is a SHAPE (MAP/M/s+G), not a presence: a patience
% declared on a two-queue model is refused by solver_mam_retrial, so the same
% predicate decides here (MAM_RENEGING_APPLICABLE).
[isReneging, renInfo] = mam_reneging_applicable(sn);
if isRetrial || isReneging
    tf = true;
    reason = '';
    return
end

tf = false;
% qsys_is_retrial reports WHICH requirement the model missed (open model,
% single class, a bufferless station, a retrial drop rule); carrying it through
% is the difference between "this model has neither" and a usable answer. A
% model that DOES declare a patience gets the reneging analyzer's own reason.
detail = '';
if mam_has_reneging_patience(sn) && isstruct(renInfo) && isfield(renInfo,'errorMsg') ...
        && ~isempty(renInfo.errorMsg)
    detail = [' ' renInfo.errorMsg];
elseif isstruct(retInfo) && isfield(retInfo,'errorMsg') && ~isempty(retInfo.errorMsg)
    detail = [' ' retInfo.errorMsg];
end
reason = sprintf(['The ''retrial'' method needs a BMAP/PH/N/N bufferless retrial ' ...
    'topology or a MAP/M/s+G reneging shape (one Source, one FCFS Queue with exponential ' ...
    'service, one open class, a patience distribution); this model has neither.%s ' ...
    'Use options.method=''default''.'], detail);
end
