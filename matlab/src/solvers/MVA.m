classdef MVA < SolverMVA
    % MVA - Alias for SolverMVA (Mean Value Analysis solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = MVA(model, varargin)
            % MVA(MODEL, VARARGIN)
            self@SolverMVA(model, varargin{:});
        end
    end
end
