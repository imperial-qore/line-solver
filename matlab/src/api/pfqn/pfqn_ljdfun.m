function r = pfqn_ljdfun(nvec, ljdscaling, ljdcutoffs, nservers)
% R = PFQN_LJDFUN(NVEC, LJDSCALING, LJDCUTOFFS, NSERVERS)
%
% Joint-dependence scaling function for AMVA-QD
%
% nvec: per-class population vector [n1, n2, ..., nR]
% ljdscaling: cell array of linearized scaling tables per station
% ljdcutoffs: M x K matrix of cutoffs
% nservers: (optional) number of servers per station
%
% Returns: r - scaling factor vector (M x 1)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = size(ljdcutoffs, 1);
K = length(nvec);
r = ones(M, 1);
alpha = 20; % softmin parameter

for i = 1:M
    if ~isempty(ljdscaling{i})
        % Clamp population to cutoffs
        nClamped = min(nvec, ljdcutoffs(i, :));

        % Compute linearized index (1-indexed for MATLAB)
        idx = ljd_linearize(nClamped, ljdcutoffs(i, :));

        % Look up scaling value
        if idx <= length(ljdscaling{i})
            r(i) = ljdscaling{i}(idx);
        end
    end

    % Handle servers (similar to pfqn_lldfun)
    if nargin >= 4
        if isinf(nservers(i)) % delay server
            % handled in the main code differently
        else
            total_n = sum(nvec);
            r(i) = r(i) / softmin(total_n, nservers(i), alpha);
            if isnan(r(i)) % if numerical problems in soft-min
                r(i) = 1 / min(total_n, nservers(i));
            end
        end
    end
end
end
