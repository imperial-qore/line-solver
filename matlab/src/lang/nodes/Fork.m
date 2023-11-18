classdef Fork < Node
    % A node to fork jobs into siblings tasks
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        schedStrategy;
        cap;
    end
    
    methods
        %Constructor
        function self = Fork(model, name)
            % SELF = FORK(MODEL, NAME)
            
            self@Node(name);
            if(model ~= 0)
                classes = model.classes;
                self.cap = Inf;
                self.input = Buffer(classes);
                self.schedStrategy = SchedStrategy.FORK;
                self.server = ServiceTunnel();
                self.output = Forker(classes);
                self.setModel(model);
                addNode(model, self);
            end
        end
        
        function setTasksPerLink(self, nTasks)
            % SETTASKSPERLINK(NTASKS)
            
            self.output.tasksPerLink = nTasks;
        end
        
        function summary(self)
            % SUMMARY()

            line_printf('\nNode: <strong>%s</strong>',self.getName);            
            for r=1:length(self.output.outputStrategy)
                line_printf('Routing %s: %s',self.model.classes{r}.name,self.output.outputStrategy{r}{2});
            end
        end
    end
    
end
