function bool = mmdp_isfeasible(Q, R)
% MMDP_ISFEASIBLE Check if (Q, R) defines a valid MMDP
%
% bool = MMDP_ISFEASIBLE(Q, R)
%
% Checks if the given matrices define a valid Markov-Modulated
% Deterministic Process (MMDP).
%
% Requirements:
% - Q must be a valid generator (square, row sums = 0, proper signs)
% - R must be diagonal with non-negative entries
% - Q and R must have compatible dimensions
%
% Parameters:
%   Q (matrix): n×n generator matrix
%   R (matrix): n×n diagonal rate matrix
%
% Returns:
%   bool (logical): True if (Q, R) defines a valid MMDP
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

    tol = 1e-10;

    % Q must be square
    [n, m] = size(Q);
    if n ~= m
        bool = false;
        return;
    end

    % Q must be a valid generator
    for i = 1:n
        % Diagonal elements must be non-positive
        if Q(i,i) > tol
            bool = false;
            return;
        end
        % Off-diagonal elements must be non-negative
        for j = 1:n
            if i ~= j && Q(i,j) < -tol
                bool = false;
                return;
            end
        end
        % Row sums must be zero (or close to zero)
        if abs(sum(Q(i,:))) > tol
            bool = false;
            return;
        end
    end

    % R must be n×n
    [nr, mr] = size(R);
    if nr ~= n || mr ~= n
        bool = false;
        return;
    end

    % R must be diagonal with non-negative entries
    for i = 1:n
        % Diagonal entries must be non-negative
        if R(i,i) < -tol
            bool = false;
            return;
        end
        % Off-diagonal entries must be zero
        for j = 1:n
            if i ~= j && abs(R(i,j)) > tol
                bool = false;
                return;
            end
        end
    end

    bool = true;
end
