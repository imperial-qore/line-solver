function [APH, APHS] = aph2_fit(M1, M2, M3)
% Fits a second order phase-type distributions, given the first three
% moments. Approximate fitting is performed if the given moments are
% infeasible.
% Input:
% - M1,M2,M3: the first three moments
% Output:
% - APH: the fitted phase-type distributions
% - APHS: all feasible forms of the fitted distribution


% perform exact fitting
APHS = aph2_fitall(M1, M2, M3);

% if no solutions is found, perform approximate fitting
if isempty(APHS)
    % find feasible characteristics
    [M2a, M3a] = aph2_adjust(M1, M2, M3);
    % fit (should found at least one solution)
    APHS = aph2_fitall(M1, M2a, M3a);
    if isempty(APHS)
        error('Fitting APH(2): feasibility could not be restored');
    else
%        fprintf('Fitting APH(2): %d approximate solutions\n', length(APHS));
%        fprintf('Fitting APH(2): M2 = %f -> %f\n', M2, M2a);
%        fprintf('Fitting APH(2): M3 = %f -> %f\n', M3, M3a);
    end
else
%    fprintf('Fitting APH(2): %d exact solutions\n', length(APHS));
end

APH = APHS{1};

end