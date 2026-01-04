classdef LogTunnel < ServiceTunnel
    % A service tunnel for Logger nodes
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = LogTunnel()
            % SELF = LOGTUNNEL()
            
            self@ServiceTunnel('LogTunnel');
            self.numberOfServers = 1;
            self.serviceProcess = {};
        end
    end
    
end

