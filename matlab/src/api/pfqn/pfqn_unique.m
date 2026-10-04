%{
%{
 % @file pfqn_unique.m
 % @brief Consolidate replicated stations into unique stations with multiplicity.
%}
%}

%{
%{
 % @brief Consolidate replicated stations into unique stations with multiplicity.
 % @fn pfqn_unique(L, mu, gamma)
 % @param L Service demand matrix (M x R).
 % @param mu Load-dependent rate matrix (M x Ntot), optional.
 % @param gamma Class-dependent service rate matrix (M x R), optional.
 % @return L_unique Reduced demand matrix (M' x R) with M' <= M unique stations.
 % @return mu_unique Reduced load-dependent rates (M' x Ntot), empty if mu was empty.
 % @return gamma_unique Reduced class-dependent rates (M' x R), empty if gamma was empty.
 % @return mi Multiplicity vector (1 x M'), mi(j) = count of stations mapping to unique station j.
 % @return mapping Mapping vector (1 x M), mapping(i) = unique station index for original station i.
%}
%}
function [L_unique, mu_unique, gamma_unique, mi, mapping] = pfqn_unique(L, mu, gamma)
% PFQN_UNIQUE Consolidate replicated stations into unique stations with multiplicity
%
% [L_UNIQUE, MU_UNIQUE, GAMMA_UNIQUE, MI, MAPPING] = PFQN_UNIQUE(L, MU, GAMMA)
%
% Identifies stations with identical demand rows L(i,:) and (if present)
% identical load-dependent rates mu(i,:) or class-dependent rates gamma(i,:).
% Returns reduced matrices with only unique stations plus a multiplicity vector.
%
% Input:
%   L     - M x R demand matrix
%   mu    - M x Ntot load-dependent rate matrix (optional, pass [] if not used)
%   gamma - M x R class-dependent service rate matrix (optional, pass [] if not used)
%
% Output:
%   L_unique     - M' x R demand matrix with M' <= M unique stations
%   mu_unique    - M' x Ntot reduced load-dependent rates (empty if mu was empty)
%   gamma_unique - M' x R reduced class-dependent rates (empty if gamma was empty)
%   mi           - 1 x M' multiplicity vector
%   mapping      - 1 x M vector mapping original station i to unique station index

[M, R] = size(L);
tol = GlobalConstants.Zero();
if isempty(tol)
    tol = 1e-14;  % Default value matching ZERO_THRESHOLD in lineStart.m
end

% Handle empty optional parameters
if nargin < 2
    mu = [];
end
if nargin < 3
    gamma = [];
end

% Build fingerprint matrix for comparison
fingerprint = L;
if ~isempty(mu)
    fingerprint = [fingerprint, mu];
end
if ~isempty(gamma)
    fingerprint = [fingerprint, gamma];
end

% Find unique rows using tolerance-based comparison
mapping = zeros(1, M);
unique_idx = [];    % indices of unique stations
mi_list = [];       % multiplicities

for i = 1:M
    if mapping(i) == 0  % not yet assigned
        % This is a new unique station
        unique_idx(end+1) = i; %#ok<AGROW>
        group_idx = length(unique_idx);
        mapping(i) = group_idx;
        count = 1;

        % Find all stations identical to this one
        for j = (i+1):M
            if mapping(j) == 0  % not yet assigned
                if norm(fingerprint(i,:) - fingerprint(j,:), inf) < tol
                    mapping(j) = group_idx;
                    count = count + 1;
                end
            end
        end
        mi_list(end+1) = count; %#ok<AGROW>
    end
end

% Extract unique rows
mi = mi_list;
L_unique = L(unique_idx, :);
if ~isempty(mu)
    mu_unique = mu(unique_idx, :);
else
    mu_unique = [];
end
if ~isempty(gamma)
    gamma_unique = gamma(unique_idx, :);
else
    gamma_unique = [];
end
end
