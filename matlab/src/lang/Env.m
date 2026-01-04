classdef Env < Environment
    % ENV is a deprecated alias for Environment.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = Env(varargin)
            % SELF = ENV(VARARGIN)
            self@Environment(varargin{:});
            warning('LINE:Deprecated', 'The Env class is deprecated. Please use Environment instead.');
        end
    end
end
