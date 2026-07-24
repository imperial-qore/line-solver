classdef InfiniteServer < ServiceSection
    % A service section to delay jobs without contention
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = InfiniteServer(customerClasses)
            % SELF = INFINITESERVER(CUSTOMERCLASSES)
            
            self@ServiceSection('InfiniteServer');
            self.numberOfServers = Inf;
            self.serviceProcess = {};
            initServers(self, customerClasses);
        end
    end
    
    methods (Access = 'private')
        function initServers(self, customerClasses)
            % INITSERVERS(CUSTOMERCLASSES)
            for i = 1 : length(customerClasses)
                self.serviceProcess{1, i} = {customerClasses{i}.name, ServiceStrategy.LI, Exp(0.0)};
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
            for i=1:length(self.serviceProcess)
                if ishandle(self.serviceProcess{i}{3})
                    clone.serviceProcess{i}{3} = self.serviceProcess{i}{3}.copy;
                end
            end
        end
    end
    
end

