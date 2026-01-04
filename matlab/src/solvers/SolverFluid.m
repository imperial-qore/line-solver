classdef SolverFluid < SolverFLD
    % SolverFLD - Alias for SolverFLD (Fluid/Mean-Field Approximation solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverFluid(model, varargin)
            % SolverFluid(MODEL, VARARGIN)
            self@SolverFLD(model, varargin{:});
        end
    end
end
