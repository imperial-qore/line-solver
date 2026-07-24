classdef SSA < SolverSSA
    % SSA - Alias for SolverSSA (Stochastic State-space Analysis solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SSA(model, varargin)
            % SSA(MODEL, VARARGIN)
            self@SolverSSA(model, varargin{:});
        end
    end
end
