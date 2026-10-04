classdef Timing < ServiceSection
    % A service tunnel for general nodes
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = Timing(name)
            % SELF = SERVICETUNNEL(NAME)
            
            if nargin<1
                name = 'Timing';
            end
            self@ServiceSection(name);
            self.numberOfServers = 1;
            self.serviceProcess = {};
        end
    end
    
end

