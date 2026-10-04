classdef AUTO < SolverAuto
    % AUTO - Alias for SolverAuto (Automatic solver selection)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = AUTO(model, varargin)
            % AUTO(MODEL, VARARGIN)
            self@SolverAuto(model, varargin{:});
        end
    end
end
