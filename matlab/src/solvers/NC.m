classdef NC < SolverNC
    % NC - Alias for SolverNC (Normalizing Constant solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = NC(model, varargin)
            % NC(MODEL, VARARGIN)
            self@SolverNC(model, varargin{:});
        end
    end
end
