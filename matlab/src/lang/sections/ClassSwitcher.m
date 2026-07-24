classdef ClassSwitcher < ServiceSection
    % An abstract classs for a service section that switches job class
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        csFun;
        classes;
    end
    
    methods
        %Constructor
        function self = ClassSwitcher(classes, name)
            % SELF = CLASSSWITCHER(CLASSES, NAME)
            
            self@ServiceSection(name);
            self.classes = classes;
            self.numberOfServers = 1;
            self.serviceProcess = {};
        end
    end
    
end
