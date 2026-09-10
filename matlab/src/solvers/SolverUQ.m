classdef SolverUQ < UQ
    % SolverUQ - Alias for UQ (uncertainty quantification solver wrapper)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverUQ(model, solverFactory, varargin)
            % SOLVERUQ(MODEL, SOLVERFACTORY, VARARGIN)
            self@UQ(model, solverFactory, varargin{:});
        end
    end
end
