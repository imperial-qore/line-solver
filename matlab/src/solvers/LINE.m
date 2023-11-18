classdef LINE < SolverAuto
    % Alias for SolverAuto
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = LINE(model, varargin)
            % SELF = LINE(MODEL, VARARGIN)
            self@SolverAuto(model, varargin{:});
        end
    end

end
