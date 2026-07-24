classdef Environment < Ensemble
    % An environment model defined by a collection of network sub-models
    % coupled with an environment transition rule that selects the active
    % sub-model.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        env;
        envGraph;
        num_stages; % Expected number of stages (optional, for API consistency)
        proc; % Markovian representation of each stage transition
        holdTime; % holding times
        probEnv; % steady-stage probability of the environment
        probOrig; % probability that a request originated from phase
        resetFun; % function implementing the reset policy for queue lengths
        resetEnvRatesFun; % function implementing the reset policy for environment rates
        resetStateFun; % function implementing the reset policy for the full state-probability vector (SolverENV statevec analyzer)
        nodeFailures; % cell array of node breakdown/repair descriptors recorded by addNodeBreakdown/addNodeRepair
    end
    
    methods
        function self = Environment(name, num_stages)
            % SELF = ENVIRONMENT(NAME, NUM_STAGES)
            % NAME - Name of the environment
            % NUM_STAGES - (Optional) Expected number of stages. If provided, used for validation only.
            %              Stages can still be added/removed dynamically.
            if nargin < 2
                num_stages = 0;
            end
            self@Ensemble({});
            self.name = name;
            self.num_stages = num_stages;
            self.nodeFailures = {};
            self.envGraph = digraph();
            self.envGraph.Nodes.Model = cell(0);
            self.envGraph.Nodes.Type = cell(0);
        end
        
        function name = addStage(self, name, type, model)            
            wcfg = warning; % store warning configuration
            warning('off','MATLAB:table:RowsAddedExistingVars');
            self.envGraph = self.envGraph.addnode(name);
            warning(wcfg); % restore warning configuration
            self.envGraph.Nodes.Model{end} = model;
            self.envGraph.Nodes.Type{end} = type;
            E = height(self.envGraph.Nodes);
            if E>1
                if self.envGraph.Nodes.Model{1}.getNumberOfStatefulNodes ~= model.getNumberOfStatefulNodes
                    line_error(mfilename,'Unsupported feature. Random environment stages must map to networks with identical number of stateful nodes.');
                end
            end
            for e=E
                for h=1:E
                    self.env{e,h} = Disabled.getInstance();
                    self.env{h,e} = Disabled.getInstance();
                end
            end
            self.ensemble{E} = model;
        end
        
        function self = addTransition(self, fromName, toName, distrib, resetFun, resetEnvRatesFun, resetStateFun)
           self.envGraph = self.envGraph.addedge(fromName, toName);
           e = self.envGraph.findnode(fromName);
           h = self.envGraph.findnode(toName);
           self.env{e,h} = distrib;
           if nargin<5 %~exist('resetFun','var')
                self.resetFun{e,h} = @(q) q;
           else
                self.resetFun{e,h} = resetFun;
           end
           if nargin<6 %~exist('resetEnvRatesFun','var')
                self.resetEnvRatesFun{e,h} = @(originalDist, QExit, UExit, TExit) originalDist;
           else
                self.resetEnvRatesFun{e,h} = resetEnvRatesFun;
           end
           if nargin<7 %~exist('resetStateFun','var')
                self.resetStateFun{e,h} = @(pi) pi; % identity: carry the exit distribution unchanged
           else
                self.resetStateFun{e,h} = resetStateFun;
           end
        end
        
        function init(self)
            E = height(self.envGraph.Nodes);
            T = height(self.envGraph.Edges);
            Pemb = zeros(E);

            for t=1:T
                e = self.envGraph.findnode(self.envGraph.Edges.EndNodes{t,1});
                h = self.envGraph.findnode(self.envGraph.Edges.EndNodes{t,2});
                self.envGraph.Edges.Distribution{t} = self.env{e,h};
            end
            
            % analyse holding times
            emmap = cell(E,1);
            for e=1:E
                emmap{e} = cell(1,E);
            end
            self.holdTime = {};
            for e=1:E
                for h=1:E
                    if isa(self.env{e,h},'Disabled')
                        emmap{e}{h} = {0,0}; % multiclass MMAP representation
                    else
                        emmap{e}{h} = self.env{e,h}.getProcess; % multiclass MMAP representation
                    end
                    for j = 1:E
                        if j == h
                            emmap{e}{h}{2+j} = emmap{e}{h}{2};
                        else
                            emmap{e}{h}{2+j} = 0 * emmap{e}{h}{2};
                        end
                    end
                end
                self.holdTime{e} = emmap{e}{e};
                for h=setdiff(1:E,e)
                    self.holdTime{e}{1} = krons(self.holdTime{e}{1},emmap{e}{h}{1});
                    for j = 2:(E+2)
                        self.holdTime{e}{j} = krons(self.holdTime{e}{j},emmap{e}{h}{j});
                        completion_rates = self.holdTime{e}{j}*ones(length(self.holdTime{e}{j}),1);
                        self.holdTime{e}{j} = 0*self.holdTime{e}{j};
                        self.holdTime{e}{j}(:,1) = completion_rates;
                    end
                    self.holdTime{e} = mmap_normalize(self.holdTime{e});
                end
                count_lambda = mmap_count_lambda(self.holdTime{e}); % completiom rates for the different transitions
                Pemb(e,:) = count_lambda/sum(count_lambda);
            end
            self.proc = emmap;
            
            %
            lambda = zeros(1,E);
            A = zeros(E); I=eye(E);
            for e=1:E
                lambda(e) = 1/map_mean(self.holdTime{e});
                for h=1:E
                    A(e,h) = -lambda(e)*(I(e,h)-Pemb(e,h));
                end
            end
            
            if all(lambda>0)
                penv = ctmc_solve_reducible(A);
                self.probEnv = penv;
                self.probOrig = zeros(E);
                for e = 1:E
                    for h = 1:E
                        self.probOrig(h,e) = penv(h) * lambda(h) * Pemb(h,e);
                    end
                    if penv(e) > 0
                        self.probOrig(:,e) = self.probOrig(:,e) / sum(self.probOrig(:,e));
                    end
                end
            else
                %T = A(lambda>0, lambda>0), % subgenerator
                %t = A(lambda>0, lambda==0), % exit vector
            end
        end
        
        function env = getEnv(self)
            env = self.env;
        end
        
        function self = setEnv(self, env)
            self.envGraph = env;
        end
        
        function self = setStageName(self, stageId, name)
            self.stageNames{stageId} = name;
        end
                
        function self = setStageType(self, stageId, stageCategory)                      
            if ischar(stageCategory)
                self.stageTypes(stageId) = categorical(stageCategory);
            elseif iscategorical(stageCategory)
                self.stageTypes(stageId) = stageCategory;
            else
                line_error(mfilename,'Stage type must be of type categorical, e.g., categorical("My Semantics").');
            end
        end
        
        function ET = getStageTable(self)
            E = height(self.envGraph.Nodes);
            Stage = [];
            HoldT = {};
            type = categorical([]);
            if isempty(self.probEnv)
                self.init;
            end
            for e=1:E
                Stage(e,1) = e;
                type(e,1) = self.envGraph.Nodes.Type{e};
            end
            Prob = self.probEnv(:);
            Name = categorical(self.envGraph.Nodes.Name(:));
            Model = self.ensemble(:);
            for e=1:E
                HoldT{e,1} = self.holdTime{e};
            end
            Type = type;
            ET = Table(Stage, Name, Type, Prob, HoldT, Model);
        end

        function ET = getStageT(self)
            % GETSTAGET Short alias for getStageTable
            ET = self.getStageTable();
        end

        function RT = getReliabilityTable(self)
            % RT = GETRELIABILITYTABLE()
            % Compute system-wide reliability metrics (MTTF, MTTR, MTBF, Availability)
            %
            % Returns:
            %   RT - Table with columns: Metric, Value, Unit, Description
            %
            % Example:
            %   env = Environment('ServerEnv');
            %   env.addNodeFailureRepair(model, 'Server', Exp(0.1), Exp(1.0), Exp(0.5));
            %   env.init();
            %   reliabilityTable = env.getReliabilityTable();

            % Step 1: Initialize and validate
            if isempty(self.probEnv)
                self.init();
            end

            E = height(self.envGraph.Nodes);
            if E == 0
                line_error(mfilename, 'Environment has no stages. Add stages before computing reliability metrics.');
            end

            % Step 2: Identify stage types
            stageNames = self.envGraph.Nodes.Name;
            upIdx = find(strcmp(stageNames, 'UP'));

            if isempty(upIdx)
                line_error(mfilename, 'No UP stage found. Use addNodeBreakdown/addNodeRepair to configure breakdown/repair transitions.');
            end

            % Find all DOWN stages
            downIdx = find(startsWith(stageNames, 'DOWN_'));

            if isempty(downIdx)
                line_error(mfilename, 'No DOWN stages found. Use addNodeBreakdown/addNodeRepair to configure breakdown/repair transitions.');
            end

            % Step 3: Extract breakdown rates (UP -> DOWN_*)
            breakdownRates = [];
            for h = downIdx'
                if ~isa(self.env{upIdx, h}, 'Disabled')
                    lambda_h = 1 / self.env{upIdx, h}.getMean();
                    breakdownRates(end+1) = lambda_h;
                end
            end

            if isempty(breakdownRates)
                line_error(mfilename, 'No breakdown transitions found (UP -> DOWN_*).');
            end

            % Total failure rate (competing risks)
            lambda_total = sum(breakdownRates);
            MTTF = 1 / lambda_total;

            % Step 4: Extract repair rates (DOWN_* -> UP)
            repairRates = [];
            downProbs = [];

            for e = downIdx'
                if ~isa(self.env{e, upIdx}, 'Disabled')
                    mu_e = 1 / self.env{e, upIdx}.getMean();
                    repairRates(end+1) = mu_e;
                    downProbs(end+1) = self.probEnv(e);
                end
            end

            if isempty(repairRates)
                line_error(mfilename, 'No repair transitions found (DOWN_* -> UP).');
            end

            % Normalize probabilities over DOWN states only
            totalDownProb = sum(downProbs);
            if totalDownProb > 0
                downProbsNorm = downProbs / totalDownProb;
                % Weighted average repair time
                MTTR = sum(downProbsNorm ./ repairRates);
            else
                % Fallback: simple average if no steady-state probability
                MTTR = mean(1 ./ repairRates);
            end

            % Step 5: Compute derived metrics
            MTBF = MTTF + MTTR;

            % Availability from steady-state probabilities
            availUp = self.probEnv(upIdx);
            availDown = sum(self.probEnv(downIdx));
            Availability = availUp / (availUp + availDown);

            % Step 6: Create output table
            Metric = categorical({'MTTF'; 'MTTR'; 'MTBF'; 'Availability'});
            Value = [MTTF; MTTR; MTBF; Availability];
            Unit = categorical({'time units'; 'time units'; 'time units'; 'probability'});
            Description = categorical({
                'Mean time to failure (UP -> DOWN)';
                'Mean time to repair (DOWN -> UP)';
                'Mean time between failures (MTTF + MTTR)';
                'Steady-state probability of UP state'
            });

            RT = Table(Metric, Value, Unit, Description);
        end

        function RT = relT(self)
            % RELT Short alias for getReliabilityTable
            RT = self.getReliabilityTable();
        end

        function RT = getRelT(self)
            % GETRELT Short alias for getReliabilityTable
            RT = self.getReliabilityTable();
        end

        function RT = relTable(self)
            % RELTABLE Short alias for getReliabilityTable
            RT = self.getReliabilityTable();
        end

        function RT = getRelTable(self)
            % GETRELTABLE Short alias for getReliabilityTable
            RT = self.getReliabilityTable();
        end

        function printStageTable(self)
            % PRINTSTAGETABLE Print a formatted table showing all stages, their properties, and transitions
            %
            % Displays stage names, types, associated networks, and transition rates.
            %
            % Example:
            %   env = Environment('MyEnv');
            %   env.addStage('UP', 'operational', model1);
            %   env.addStage('DOWN', 'failed', model2);
            %   env.addTransition('UP', 'DOWN', Exp(0.1));
            %   env.printStageTable();

            E = height(self.envGraph.Nodes);
            fprintf('Stage Table:\n');
            fprintf('============\n');

            for e = 1:E
                stageName = self.envGraph.Nodes.Name{e};
                stageType = self.envGraph.Nodes.Type{e};
                model = self.envGraph.Nodes.Model{e};

                fprintf('Stage %d: %s (Type: %s)\n', e, stageName, stageType);
                if ~isempty(model)
                    fprintf('  - Network: %s\n', model.getName());
                    fprintf('  - Nodes: %d\n', model.getNumberOfNodes());
                    fprintf('  - Classes: %d\n', model.getNumberOfClasses());
                end
            end

            % Print transitions
            T = height(self.envGraph.Edges);
            if T > 0
                fprintf('\nTransitions:\n');
                for t = 1:T
                    fromName = self.envGraph.Edges.EndNodes{t, 1};
                    toName = self.envGraph.Edges.EndNodes{t, 2};
                    e = self.envGraph.findnode(fromName);
                    h = self.envGraph.findnode(toName);
                    if ~isa(self.env{e, h}, 'Disabled')
                        rate = 1 / self.env{e, h}.getMean();
                        fprintf('  %s -> %s: rate = %.4f\n', fromName, toName, rate);
                    end
                end
            end
        end

        function self = addNodeBreakdown(self, baseModel, nodeOrName, breakdownDist, downServiceDist, varargin)
            % SELF = ADDNODEBREAKDOWN(BASEMODEL, NODEORNAME, BREAKDOWNDIST, DOWNSERVICEDIST, RESETFUN)
            % Adds UP and DOWN stages for a node that can break down and repair
            %
            % Parameters:
            %   baseModel - The base network model with normal (UP) service rates
            %   nodeOrName - Node object or name of the node that can break down
            %   breakdownDist - Distribution for time until breakdown (UP->DOWN transition)
            %   downServiceDist - Service distribution when the node is down
            %   resetFun - (Optional) Reset policy for queue lengths on breakdown.
            %              Either a function handle @(q) -> q, or one of the named
            %              policies 'keep' (identity) and 'clear' (empty the queues).
            %              Default: 'keep'. Only named policies are serializable.
            %
            % Example:
            %   model = Network('MyNetwork');
            %   queue = Queue(model, 'Server1', SchedStrategy.FCFS);
            %   class = ClosedClass(model, 'Jobs', 10, queue, 0);
            %   queue.setService(class, Exp(2)); % UP service rate
            %
            %   env = Environment('ServerEnv');
            %   env.addNodeBreakdown(model, 'Server1', Exp(0.1), Exp(0.5));
            %   % Or using node object:
            %   env.addNodeBreakdown(model, queue, Exp(0.1), Exp(0.5));

            % Extract node name if a Node object is passed
            if isa(nodeOrName, 'Node')
                nodeName = nodeOrName.name;
            else
                nodeName = nodeOrName;
            end

            if nargin < 6
                resetSpec = 'keep';
            else
                resetSpec = varargin{1};
            end
            [resetFun, resetName] = Environment.resolveResetPolicy(resetSpec);

            % Create UP stage (if this is the first call)
            if height(self.envGraph.Nodes) == 0
                upModel = baseModel.copy();
                self.addStage('UP', 'operational', upModel);
            end

            % Create DOWN stage with modified service rate for the specified node
            downModel = baseModel.copy();
            nodes = downModel.getNodes();
            nodeIdx = [];
            for i = 1:length(nodes)
                if strcmp(nodes{i}.name, nodeName)
                    nodeIdx = i;
                    break;
                end
            end

            if isempty(nodeIdx)
                line_error(mfilename, sprintf('Node "%s" not found in the base model.', nodeName));
            end

            % Update service distribution for the down node
            classes = downModel.getClasses();
            for c = 1:length(classes)
                nodes{nodeIdx}.setService(classes{c}, downServiceDist);
            end

            % Add DOWN stage
            downStageName = sprintf('DOWN_%s', nodeName);
            self.addStage(downStageName, 'failed', downModel);

            % Add breakdown transition (UP -> DOWN)
            self.addTransition('UP', downStageName, breakdownDist, resetFun);

            % Record the descriptor so that the breakdown can be serialized
            % declaratively (linemodel_save "nodeFailures"). The expanded
            % stages/transitions carry the structure, but the reset policy is a
            % function and is only recoverable from this record.
            nf = struct();
            nf.node = nodeName;
            nf.breakdown = breakdownDist;
            nf.downService = downServiceDist;
            nf.repair = [];
            nf.breakdownResetPolicy = resetName;
            nf.repairResetPolicy = '';
            self.nodeFailures{end+1} = nf;
        end

        function self = addNodeRepair(self, nodeOrName, repairDist, varargin)
            % SELF = ADDNODEREPAIR(NODEORNAME, REPAIRDIST, RESETFUN)
            % Adds repair transition from DOWN to UP stage for a previously added breakdown
            %
            % Parameters:
            %   nodeOrName - Node object or name of the node that can be repaired
            %   repairDist - Distribution for repair time (DOWN->UP transition)
            %   resetFun - (Optional) Reset policy for queue lengths on repair.
            %              Either a function handle @(q) -> q, or one of the named
            %              policies 'keep' (identity) and 'clear' (empty the queues).
            %              Default: 'keep'. Only named policies are serializable.
            %
            % Example:
            %   env.addNodeRepair('Server1', Exp(1.0));
            %   % Or using node object:
            %   env.addNodeRepair(queue, Exp(1.0));

            % Extract node name if a Node object is passed
            if isa(nodeOrName, 'Node')
                nodeName = nodeOrName.name;
            else
                nodeName = nodeOrName;
            end

            if nargin < 4
                resetSpec = 'keep';
            else
                resetSpec = varargin{1};
            end
            [resetFun, resetName] = Environment.resolveResetPolicy(resetSpec);

            downStageName = sprintf('DOWN_%s', nodeName);

            % Verify DOWN stage exists
            if isempty(self.envGraph.findnode(downStageName))
                line_error(mfilename, sprintf('DOWN stage for node "%s" not found. Call addNodeBreakdown first.', nodeName));
            end

            % Add repair transition (DOWN -> UP)
            self.addTransition(downStageName, 'UP', repairDist, resetFun);

            % Complete the descriptor recorded by addNodeBreakdown.
            idx = self.findNodeFailure(nodeName);
            if idx > 0
                self.nodeFailures{idx}.repair = repairDist;
                self.nodeFailures{idx}.repairResetPolicy = resetName;
            end
        end

        function self = registerNodeFailure(self, nodeName, breakdownDist, repairDist, downServiceDist, breakdownPolicy, repairPolicy)
            % SELF = REGISTERNODEFAILURE(NODENAME, BREAKDOWNDIST, REPAIRDIST, DOWNSERVICEDIST, BREAKDOWNPOLICY, REPAIRPOLICY)
            % Attach a node breakdown/repair descriptor to stages that already exist.
            %
            % This is the counterpart of addNodeBreakdown/addNodeRepair for the case
            % where the UP and DOWN_<node> stages and their transitions have already
            % been built (for instance by linemodel_load reading the expanded
            % stages/transitions form). It records the descriptor and applies the
            % queue-length reset policies, which the expanded form cannot carry.
            downStageName = sprintf('DOWN_%s', nodeName);
            upIdx = self.envGraph.findnode('UP');
            downIdx = self.envGraph.findnode(downStageName);
            if isempty(upIdx) || upIdx == 0
                line_error(mfilename, sprintf(['Cannot register a node failure on "%s": no UP stage is ' ...
                    'defined in this environment.'], nodeName));
            end
            if isempty(downIdx) || downIdx == 0
                line_error(mfilename, sprintf(['Cannot register a node failure on "%s": no "%s" stage is ' ...
                    'defined in this environment.'], nodeName, downStageName));
            end

            [breakdownFun, breakdownName] = Environment.resolveResetPolicy(breakdownPolicy);
            self.resetFun{upIdx, downIdx} = breakdownFun;
            repairName = '';
            if ~isempty(repairDist)
                [repairFun, repairName] = Environment.resolveResetPolicy(repairPolicy);
                self.resetFun{downIdx, upIdx} = repairFun;
            end

            nf = struct();
            nf.node = nodeName;
            nf.breakdown = breakdownDist;
            nf.downService = downServiceDist;
            nf.repair = repairDist;
            nf.breakdownResetPolicy = breakdownName;
            nf.repairResetPolicy = repairName;
            idx = self.findNodeFailure(nodeName);
            if idx > 0
                self.nodeFailures{idx} = nf;
            else
                self.nodeFailures{end+1} = nf;
            end
        end

        function idx = findNodeFailure(self, nodeName)
            % IDX = FINDNODEFAILURE(NODENAME)
            % Index of the node-failure descriptor for NODENAME, or 0 if absent.
            idx = 0;
            for i = 1:length(self.nodeFailures)
                if strcmp(self.nodeFailures{i}.node, nodeName)
                    idx = i;
                    return;
                end
            end
        end

        function self = addNodeFailureRepair(self, baseModel, nodeOrName, breakdownDist, repairDist, downServiceDist, varargin)
            % SELF = ADDNODEFAILUREREPAIR(BASEMODEL, NODEORNAME, BREAKDOWNDIST, REPAIRDIST, DOWNSERVICEDIST, RESETBREAKDOWN, RESETREPAIR)
            % Convenience method to add both breakdown and repair for a node
            %
            % Parameters:
            %   baseModel - The base network model with normal (UP) service rates
            %   nodeOrName - Node object or name of the node that can break down and repair
            %   breakdownDist - Distribution for time until breakdown
            %   repairDist - Distribution for repair time
            %   downServiceDist - Service distribution when the node is down
            %   resetBreakdown - (Optional) Reset function for breakdown transition
            %   resetRepair - (Optional) Reset function for repair transition
            %
            % Example:
            %   env = Environment('ServerEnv');
            %   env.addNodeFailureRepair(model, 'Server1', Exp(0.1), Exp(1.0), Exp(0.5));
            %   % Or using node object:
            %   env.addNodeFailureRepair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5));

            % Extract node name if a Node object is passed
            if isa(nodeOrName, 'Node')
                nodeName = nodeOrName.name;
            else
                nodeName = nodeOrName;
            end

            resetBreakdown = 'keep';
            resetRepair = 'keep';

            if nargin >= 7
                resetBreakdown = varargin{1};
            end
            if nargin >= 8
                resetRepair = varargin{2};
            end

            self.addNodeBreakdown(baseModel, nodeName, breakdownDist, downServiceDist, resetBreakdown);
            self.addNodeRepair(nodeName, repairDist, resetRepair);
        end

        function self = setBreakdownResetPolicy(self, nodeOrName, resetFun)
            % SELF = SETBREAKDOWNRESETPOLICY(NODEORNAME, RESETFUN)
            % Update the reset policy for breakdown transitions (UP -> DOWN) of a node
            %
            % Parameters:
            %   nodeOrName - Node object or name of the node
            %   resetFun - Reset policy for queue lengths on breakdown. Either a
            %              function handle, or one of the named policies 'keep'
            %              (identity) and 'clear' (empty the queues). Only named
            %              policies are serializable.
            %              Example: @(q) 0*q to clear queues, @(q) q to keep jobs
            %
            % Example:
            %   env.setBreakdownResetPolicy('Server1', @(q) 0*q);
            %   % Or using node object:
            %   env.setBreakdownResetPolicy(queue, 'clear');

            % Extract node name if a Node object is passed
            if isa(nodeOrName, 'Node')
                nodeName = nodeOrName.name;
            else
                nodeName = nodeOrName;
            end

            [resetFun, resetName] = Environment.resolveResetPolicy(resetFun);
            downStageName = sprintf('DOWN_%s', nodeName);

            % Find UP and DOWN stage indices
            upIdx = self.envGraph.findnode('UP');
            downIdx = self.envGraph.findnode(downStageName);

            if isempty(upIdx) || upIdx == 0
                line_error(mfilename, 'UP stage not found. Call addNodeBreakdown first.');
            end
            if isempty(downIdx) || downIdx == 0
                line_error(mfilename, sprintf('DOWN stage for node "%s" not found. Call addNodeBreakdown first.', nodeName));
            end

            % Update the reset function for the breakdown transition (UP -> DOWN)
            self.resetFun{upIdx, downIdx} = resetFun;
            idx = self.findNodeFailure(nodeName);
            if idx > 0
                self.nodeFailures{idx}.breakdownResetPolicy = resetName;
            end
        end

        function self = setRepairResetPolicy(self, nodeOrName, resetFun)
            % SELF = SETREPAIRRESETPOLICY(NODEORNAME, RESETFUN)
            % Update the reset policy for repair transitions (DOWN -> UP) of a node
            %
            % Parameters:
            %   nodeOrName - Node object or name of the node
            %   resetFun - Reset policy for queue lengths on repair. Either a
            %              function handle, or one of the named policies 'keep'
            %              (identity) and 'clear' (empty the queues). Only named
            %              policies are serializable.
            %              Example: @(q) q to keep jobs, @(q) 0*q to clear queues
            %
            % Example:
            %   env.setRepairResetPolicy('Server1', @(q) q);
            %   % Or using node object:
            %   env.setRepairResetPolicy(queue, 'keep');

            % Extract node name if a Node object is passed
            if isa(nodeOrName, 'Node')
                nodeName = nodeOrName.name;
            else
                nodeName = nodeOrName;
            end

            [resetFun, resetName] = Environment.resolveResetPolicy(resetFun);
            downStageName = sprintf('DOWN_%s', nodeName);

            % Find UP and DOWN stage indices
            upIdx = self.envGraph.findnode('UP');
            downIdx = self.envGraph.findnode(downStageName);

            if isempty(upIdx) || upIdx == 0
                line_error(mfilename, 'UP stage not found. Call addNodeBreakdown first.');
            end
            if isempty(downIdx) || downIdx == 0
                line_error(mfilename, sprintf('DOWN stage for node "%s" not found. Call addNodeBreakdown first.', nodeName));
            end

            % Update the reset function for the repair transition (DOWN -> UP)
            self.resetFun{downIdx, upIdx} = resetFun;
            idx = self.findNodeFailure(nodeName);
            if idx > 0
                self.nodeFailures{idx}.repairResetPolicy = resetName;
            end
        end
    end

    methods (Static)
        function [resetFun, resetName] = resolveResetPolicy(spec)
            % [RESETFUN, RESETNAME] = RESOLVERESETPOLICY(SPEC)
            % Resolve a queue-length reset policy given either a named policy or a
            % function handle.
            %
            % Named policies (the only serializable ones):
            %   'keep'  - carry the queue lengths across the transition, @(q) q
            %   'clear' - empty the queues on the transition, @(q) 0*q
            %
            % A function handle is returned unchanged and reported as 'custom':
            % an arbitrary reset function cannot be reproduced from JSON.
            if ischar(spec) || isstring(spec)
                switch lower(char(spec))
                    case 'keep'
                        resetFun = @(q) q;
                        resetName = 'keep';
                    case 'clear'
                        resetFun = @(q) 0*q;
                        resetName = 'clear';
                    otherwise
                        line_error(mfilename, sprintf(['Unknown reset policy "%s". Use ''keep'', ''clear'', ' ...
                            'or a function handle @(q) -> q.'], char(spec)));
                end
            elseif isa(spec, 'function_handle')
                resetFun = spec;
                resetName = 'custom';
            else
                line_error(mfilename, ['Reset policy must be a function handle @(q) -> q, or one of the ' ...
                    'named policies ''keep'' and ''clear''.']);
            end
        end
    end
end
