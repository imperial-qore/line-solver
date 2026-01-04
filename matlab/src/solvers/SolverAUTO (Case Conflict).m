classdef SolverAUTO < SolverAuto
    % SolverAUTO - Alias for SolverAuto (Automatic solver selection)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverAUTO(model, varargin)
            % SolverAUTO(MODEL, VARARGIN)
            self@SolverAuto(model, varargin{:});
        end
    end
end
