classdef LN < SolverLN
    % LN - Alias for SolverLN (Layered Network solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = LN(model, varargin)
            % LN(MODEL, VARARGIN)
            self@SolverLN(model, varargin{:});
        end
    end
end
