classdef PollingServer < ServiceSection
    % A service section for serving input buffers cyclically
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = PollingServer(classes)
            % SELF = SERVER(CLASSES)
            
            self@ServiceSection('PollingServer');
            self.numberOfServers = 1;
            self.serviceProcess = {};
            initServers(self, classes);
        end
    end
    
    methods (Access = 'private')
        function initServers(self, classes)
            % INITSERVERS(CLASSES)
            
            for i = 1 : length(classes)
                self.serviceProcess{1, i} = {classes{i}.name, ServiceStrategy.LI, Exp(0.0)};
            end
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            for i = 1 : length(self.serviceProcess)
                if ishandle(clone.serviceProcess{1,i}{3})
                    % this is a problem if one modifies the classes in the
                    % model because the one below is not an handle so it
                    % will not be modified
                    clone.serviceProcess{1,i}{3} = self.serviceProcess{1,i}{3}.copy;
                end
            end
        end
    end
    
end

