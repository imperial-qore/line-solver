classdef SharedServer < ServiceSection
    % A service section for serving jobs with processor sharing type
    % policies
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        %Constructor
        function self = SharedServer(customerClasses)
            % SELF = SHAREDSERVER(CUSTOMERCLASSES)
            
            self@ServiceSection('SharedServer');
            self.numberOfServers = 1;
            self.serviceProcess = {};
            initServers(self, customerClasses); %[JobClass(), ServiceStrategy.LI, Distribution('exp')];
        end
    end
    
    methods (Access = 'private')
        function initServers(self, customerClasses)
            % INITSERVERS(CUSTOMERCLASSES)
            
            for i = 1 : length(customerClasses),
                self.serviceProcess{1, i} = {customerClasses{1, i}.name, ServiceStrategy.LI, Exponential()};
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

