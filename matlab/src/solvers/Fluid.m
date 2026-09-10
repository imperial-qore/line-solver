classdef Fluid < SolverFluid
    % Fluid - Alias for SolverFluid (Fluid/Mean-Field Approximation solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = Fluid(model, varargin)
            % FLUID(MODEL, VARARGIN)
            self@SolverFluid(model, varargin{:});
        end
    end
end
