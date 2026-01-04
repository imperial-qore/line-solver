classdef DES < SolverDES
    % DES - Alias for SolverDES (Discrete Event Simulation solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = DES(model, varargin)
            % DES(MODEL, VARARGIN)
            self@SolverDES(model, varargin{:});
        end
    end
end
