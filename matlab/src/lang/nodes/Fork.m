classdef Fork < Node
    % Fork Job splitting node for parallel processing models
    %
    % Fork is a specialized node that splits incoming jobs into multiple sibling
    % tasks that can be processed in parallel by downstream nodes. Each job
    % arriving at a Fork node is replicated into multiple parallel tasks that
    % must later be synchronized using a corresponding Join node.
    %
    % @brief Job splitting node that creates parallel sibling tasks from incoming jobs
    %
    % Key characteristics:
    % - Splits each job into multiple parallel tasks
    % - Works in conjunction with Join nodes for synchronization
    % - Supports different fork strategies and task distributions
    % - Essential for modeling parallel processing systems
    % - No service delay - instantaneous job splitting
    %
    % Fork nodes are commonly used for:
    % - Parallel processing models
    % - Fork-join queueing networks
    % - Multi-threaded system modeling
    % - Task decomposition scenarios
    % - Distributed computing models
    %
    % Example:
    % @code
    % fork = Fork(model, 'TaskSplitter');
    % % Jobs entering this fork will be split into parallel tasks
    % % Must be paired with a Join node for proper synchronization
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        schedStrategy;
        cap;
    end
    
    methods
        %Constructor
        function self = Fork(model, name)
            % FORK Create a Fork node instance
            %
            % @brief Creates a Fork node for splitting jobs into parallel tasks
            % @param model Network model to add the fork to
            % @param name String identifier for the fork node
            % @return self Fork instance configured for the given model
            %
            % The constructor initializes the Fork node with appropriate buffers,
            % service tunnels, and forker output components. Fork nodes have no
            % service delay and immediately split incoming jobs into parallel tasks.
            
            self@Node(name);
            if model.isMatlabNative()
                if(model ~= 0)
                    classes = model.getClasses();
                    self.cap = Inf;
                    self.input = Buffer(classes);
                    self.schedStrategy = SchedStrategy.FORK;
                    self.server = ServiceTunnel();
                    self.output = Forker(classes);
                    self.setModel(model);
                    model.addNode(self);
                end
            elseif model.isJavaNative()
                self.setModel(model);
                self.obj = jline.lang.nodes.Fork(model.obj, name);
                self.index = model.obj.getNodeIndex(self.obj);
            end
        end
        
        function setTasksPerLink(self, nTasks)
            % SETTASKSPERLINK Configure number of tasks per output link
            %
            % Sets the number of tasks sent out on each outgoing link. By default,
            % a Fork node sends exactly one task per outgoing link. This method
            % allows configuring the Fork to send multiple identical tasks on each
            % link. The total number of tasks created will be:
            % (number of outgoing links) × tasksPerLink.
            %
            % Solver compatibility for tasksPerLink > 1:
            %   - SolverJMT: Fully supported - simulation handles multiple tasks correctly
            %   - SolverDES: Fully supported - simulation handles multiple tasks correctly
            %   - SolverMVA (H-T method): Not supported - throws error
            %   - SolverMVA (MMT method): Supported - analytical approximation
            %
            % @param nTasks Number of tasks per link (default: 1)

            if nTasks ~= 1
                line_warning(mfilename, 'The setTasksPerLink feature is experimental and results may be inaccurate for analytical solvers.');
            end
            self.output.tasksPerLink = nTasks;
        end
        
        function summary(self)
            % SUMMARY Display fork node configuration summary
            %
            % @brief Prints a summary of the fork node's routing configuration

            line_printf('\nNode: <strong>%s</strong>',self.getName);            
            for r=1:length(self.output.outputStrategy)
                classes = self.model.getClasses();
                line_printf('Routing %s: %s',classes{r}.name,self.output.outputStrategy{r}{2});
            end
        end
    end
    
end
