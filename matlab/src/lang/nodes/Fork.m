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

    properties (Hidden)
        % state-carrying slots used only on FJ tag-augmented model copies
        % (ModelAdapter.fjtag), where the Fork is treated as a stateful
        % node holding the parent job for one vanishing state; they mirror
        % the StatefulNode API without changing the class ancestry
        state = [];
        statePrior = [];
        space = {};
    end

    methods(Hidden)
        function prior = getStatePrior(self)
            % PRIOR = GETSTATEPRIOR()
            prior = self.statePrior;
        end

        function self = setStatePrior(self, prior)
            % SELF = SETSTATEPRIOR(PRIOR)
            self.statePrior = prior(:);
        end

        function self = setState(self, state)
            % SELF = SETSTATE(STATE)
            self.state = state;
        end

        function state = getState(self)
            % STATE = GETSTATE()
            state = self.state;
        end

        function state = getStateSpace(self)
            % STATE = GETSTATESPACE()
            state = self.space;
        end

        function self = setStateSpace(self, space)
            % SELF = SETSTATESPACE(SPACE)
            self.space = space;
        end

        function self = resetStateSpace(self)
            % SELF = RESETSTATESPACE()
            self.space = {};
        end
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
        
        function setTasksPerLink(self, nTasks, varargin)
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
            %   - SolverLDES: Fully supported - simulation handles multiple tasks correctly
            %   - SolverMVA (H-T method): Not supported - throws error
            %   - SolverMVA/SolverNC (MMT method): Supported - the auxiliary open
            %     class carries the load of all (links x nTasks) siblings and the
            %     join synchronises on the order statistic of that many branch
            %     times, each branch replicated nTasks times. That is the same
            %     approximation an ordinary fork-join gets, but the warning below
            %     still stands for the per-link and DISTRIBUTION forms, where the
            %     analytical solvers see only the mean fanout.
            %
            % SETTASKSPERLINK(JOBCLASS, NTASKS) sets it for one class only,
            % leaving every other class on the node-wide value.
            %
            % SETTASKSPERLINK(JOBCLASS, NTASKS, DESTNODE) sets it for the link
            % towards DESTNODE only, leaving the other links alone.
            %
            % @param nTasks Number of tasks per link (default: 1)

            if nargin > 2
                % (jobclass, nTasks [, destNode]) form
                jobclass = nTasks;
                nTasks = varargin{1};
                if length(varargin) > 1
                    destName = varargin{2}.getName();
                else
                    destName = '';
                end
                self.output.tasksPerLinkByDest(end+1) = struct('dest',destName, ...
                    'class',jobclass.index,'value',nTasks);
                if nTasks ~= 1
                    line_warning(mfilename, 'The setTasksPerLink feature is experimental and results may be inaccurate for analytical solvers.');
                end
                return
            end
            if nTasks ~= 1
                line_warning(mfilename, 'The setTasksPerLink feature is experimental and results may be inaccurate for analytical solvers.');
            end
            self.output.tasksPerLink = nTasks;
        end

        function setTasksPerLinkDistribution(self, jobclass, dist, destNode)
            % SETTASKSPERLINKDISTRIBUTION Configure a RANDOM number of tasks per link
            %
            % SETTASKSPERLINKDISTRIBUTION(JOBCLASS, DIST) makes the number of
            % tasks emitted on each outgoing link a draw from DIST, a
            % DiscreteSampler over a non-negative integer support, redrawn
            % independently for every link and every forked job. This is the
            % variable forking level of JMT's JobsPerLinkDis.
            %
            % SETTASKSPERLINKDISTRIBUTION(JOBCLASS, DIST, DESTNODE) restricts it
            % to the link towards DESTNODE, leaving the other links alone.
            %
            % Exact under SolverJMT and SolverLDES, which draw the degree at the
            % fork epoch. The analytical solvers see E[DIST]: SolverMVA's MMT
            % method uses it as the mean fanout, and SolverNC/SolverCTMC obtain
            % the visit ratios from sn_fj_visits_spn.
            %
            % @param jobclass Job class the distribution applies to
            % @param dist DiscreteSampler over the tasks-per-link support
            % @param destNode Optional destination node restricting the link

            if ~isa(dist,'DiscreteDistribution')
                line_error(mfilename, 'The tasks-per-link distribution must be a DiscreteDistribution (e.g. DiscreteSampler).');
            end
            if nargin < 4
                destName = '';
            else
                destName = destNode.getName();
            end
            self.output.tasksPerLinkDist(end+1) = struct('dest',destName, ...
                'class',jobclass.index,'dist',dist);
        end

        function setBranchProbability(self, jobclass, destNode, prob)
            % SETBRANCHPROBABILITY Activate an outgoing branch only with probability PROB
            %
            % SETBRANCHPROBABILITY(JOBCLASS, DESTNODE, PROB) makes the branch
            % towards DESTNODE fire with probability PROB for jobs of JOBCLASS,
            % and emit nothing otherwise. The branches are activated
            % independently, so the number of siblings a job produces is random
            % even when the tasks per link are deterministic.
            %
            % The matched Join must be told what to wait for: with a standard
            % join a job that skipped a branch would block forever, so a fork
            % with any branch probability below one requires JoinStrategy.PARTIAL
            % (see Join.setRequired) or a quorum.
            %
            % @param jobclass Job class the probability applies to
            % @param destNode Destination node of the branch
            % @param prob Activation probability in [0,1]

            if prob < 0 || prob > 1
                line_error(mfilename, 'A branch activation probability must lie in [0,1].');
            end
            self.output.branchProb(end+1) = struct('dest',destNode.getName(), ...
                'class',jobclass.index,'value',prob);
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
