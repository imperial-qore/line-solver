function demandEst = infer_minps(model, node, rt, class, ql)
% INFER_MINPS MINPS demand estimation method.
%
% Runs both MLPS and RPS estimators and selects the one with the
% smaller mean demand estimate.
%
% Inputs:
%   model  - LINE Network model with delay rates set and queue rates to estimate
%   node   - PS queue node (Station object)
%   rt     - response time samples (column vector, n x 1)
%   class  - class of each sample (column vector, n x 1)
%   ql     - queue lengths at arrival (n x R matrix, per-class)
%
% Returns:
%   demandEst - 1 x R vector of estimated mean service demands
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

V = node.getNumberOfServers();

% Run MLPS
demandEstMLPS = infer_mlps(model, node, rt, class, ql);
% Run RPS
demandEstRPS = infer_rps(rt, class, ql, V);

% Choose the smallest mean result
if mean(demandEstMLPS) < mean(demandEstRPS)
    demandEst = demandEstMLPS;
else
    demandEst = demandEstRPS;
end

end
