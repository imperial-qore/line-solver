classdef Entry < LayeredNetworkElement
    % An entry point of service for a Task.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        parent;
        replyActivity = {};
        openArrivalRate;
        scheduling = [];
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Entry(model, name)
            % SELF = ENTRY(MODEL, NAME)
            name = char(name);
            if nargin<2 %~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
            self@LayeredNetworkElement(name);
            self.openArrivalRate = 0.0;
            model.entries{end+1} = self;
            self.model = model;
        end
        
        function self = on(self, parent)
            % SELF = ON(SELF, PARENT)
            
            parent.addEntry(self);
            self.parent = parent;
        end
        
    end
    
end
