%{
%{
 % @file pfqn_expand.m
 % @brief Expand per-station metrics from reduced model to original dimensions.
%}
%}

%{
%{
 % @brief Expand per-station metrics from reduced model to original dimensions.
 % @fn pfqn_expand(QN, UN, CN, mapping, M_original)
 % @param QN Queue lengths from reduced model (M' x R).
 % @param UN Utilizations from reduced model (M' x R).
 % @param CN Cycle times from reduced model (M' x R).
 % @param mapping Mapping vector from pfqn_unique (1 x M), mapping(i) = unique station index.
 % @param M_original Original number of stations M.
 % @return QN_full Queue lengths in original dimensions (M x R).
 % @return UN_full Utilizations in original dimensions (M x R).
 % @return CN_full Cycle times in original dimensions (M x R).
%}
%}
function [QN_full, UN_full, CN_full] = pfqn_expand(QN, UN, CN, mapping, M_original)
% PFQN_EXPAND Expand per-station metrics from reduced model to original dimensions
%
% [QN_FULL, UN_FULL, CN_FULL] = PFQN_EXPAND(QN, UN, CN, MAPPING, M_ORIGINAL)
%
% Expands performance metrics computed on a reduced model (with unique stations)
% back to the original model dimensions by replicating values according to mapping.
%
% Input:
%   QN         - M' x R queue lengths from reduced model
%   UN         - M' x R utilizations from reduced model
%   CN         - M' x R cycle times from reduced model
%   mapping    - 1 x M vector from pfqn_unique (mapping(i) = unique station index)
%   M_original - original number of stations M
%
% Output:
%   QN_full    - M x R queue lengths in original dimensions
%   UN_full    - M x R utilizations in original dimensions
%   CN_full    - M x R cycle times in original dimensions

R = size(QN, 2);
QN_full = zeros(M_original, R);
UN_full = zeros(M_original, R);
CN_full = zeros(M_original, R);

for i = 1:M_original
    unique_idx = mapping(i);
    QN_full(i, :) = QN(unique_idx, :);
    UN_full(i, :) = UN(unique_idx, :);
    CN_full(i, :) = CN(unique_idx, :);
end
end
