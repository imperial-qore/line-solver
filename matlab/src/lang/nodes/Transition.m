classdef Transition < StatefulNode
    % A class for a stochastic Petri net transition
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        enablingConditions;
        inhibitingConditions;
        modes;
        modeNames;
        numberOfServers;
        timingStrategies;
        distributions;
        firingPriorities;
        firingWeights;
        firingOutcomes;
        firingRateDependence;
        cap;
    end

    methods
        function self = Transition(model,name)
            % TRANSITION(MODEL, NAME)

            self@StatefulNode(name);
            if model.isMatlabNative()
                classes = model.getClasses();
                self.input = Enabling(classes);
                self.output = Firing(classes);
                self.cap = Inf; % Compatible with other nodes

                self.setModel(model);
                self.model.addNode(self);

                self.server = Timing();

                self.enablingConditions = [];
                self.inhibitingConditions = [];
                self.modeNames = {};
                self.numberOfServers = [];
                self.timingStrategies = [];
                self.distributions = {};
                self.firingPriorities = [];
                self.firingWeights = [];
                self.firingOutcomes = [];
                self.firingRateDependence = {};
            elseif model.isJavaNative()
                self.setModel(model);
                self.obj = jline.lang.nodes.Transition(model.obj, name);
                self.index = model.obj.getNodeIndex(self.obj);
            end
        end

        function self = init(self)
            % SELF = INIT()

            nclasses = length(self.model.getClasses());
            nnodes = length(self.model.getNodes());

            self.enablingConditions = cell(self.getNumberOfModes);
            self.inhibitingConditions = cell(self.getNumberOfModes);
            self.firingOutcomes = cell(self.getNumberOfModes);
            for m=1:self.getNumberOfModes
                self.enablingConditions{m} = zeros(nnodes,nclasses);
                self.inhibitingConditions{m} = Inf*ones(nnodes,nclasses); % Inf = no inhibitor arc (never inhibits); zeros would inhibit at >=0 tokens, i.e. always
                self.firingOutcomes{m} = zeros(nnodes,nclasses);
            end
            self.numberOfServers = ones(1,self.getNumberOfModes);
            self.timingStrategies = repmat(TimingStrategy.TIMED,1,self.getNumberOfModes);
            self.firingWeights = ones(1,self.getNumberOfModes);
            self.firingPriorities = ones(1,self.getNumberOfModes);
            self.distributions = cell(1, self.getNumberOfModes);
            self.distributions(:) = {Exp(1)};
            % Empty == unit multiplier (marking-independent firing rate).
            self.firingRateDependence = cell(1, self.getNumberOfModes);
        end

        function mode = addMode(self, modeName)
            nclasses = length(self.model.getClasses());
            nnodes = length(self.model.getNodes());
            self.modeNames{end+1} = modeName;
            self.enablingConditions{end+1} = zeros(nnodes,nclasses);
            self.inhibitingConditions{end+1} = Inf*ones(nnodes,nclasses);
            self.numberOfServers(end+1) = 1;
            self.timingStrategies(end+1) = TimingStrategy.TIMED;
            self.firingWeights(end+1) = 1.0;
            self.firingPriorities(end+1) = 1.0;
            self.distributions{end+1} = Exp(1);
            self.firingOutcomes{end+1} = zeros(nnodes,nclasses);
            self.firingRateDependence{end+1} = [];
            mode = Mode(self,modeName);
            self.modes{end+1} = mode;
        end

        function self = setEnablingConditions(self, mode, class, inputNode, enablingCondition)
            % SELF = SETENABLINGCONDITIONS(MODE, CLASS, NODE, ENABLINGCONDITIONS)

            if isa(inputNode, 'Place')
                inputNode = self.model.getNodeIndex(inputNode.name);
                self.enablingConditions{mode}(inputNode,class) = enablingCondition;
            else
                error('Node must be a Place node.');
            end
        end

        function self = setInhibitingConditions(self, mode, class, inputNode, inhibitingCondition)
            % SELF = SETINHIBITINGCONDITIONS(MODE, CLASS, NODE, INHIBITINGCONDITIONS)

            if isa(inputNode, 'Place')
                inputNode = self.model.getNodeIndex(inputNode.name);
                self.inhibitingConditions{mode}(inputNode,class) = inhibitingCondition;
            else
                error('Node must be a Place node.');
            end
        end

        function self = setModeNames(self, mode, modeName)
            % SELF = SETMODENAMES(MODE, MODENAMES)

            self.modeNames{mode} = modeName;
        end

        function self = setNumberOfServers(self, mode, numberOfServers)
            % SELF = SETNUMBEROFSERVERS(MODE, NUMOFSERVERS)

            self.numberOfServers(mode) = numberOfServers;
        end

        function self = setTimingStrategy(self, mode, timingStrategy)
            % SELF = SETTIMINGSTRATEGY(MODE, TIMINGSTRATEGY)

            self.timingStrategies(mode) = timingStrategy;
        end

        function self = setFiringPriorities(self, mode, firingPriority)
            % SELF = SETFIRINGPRIORITIES(MODE, FIRINGPRIORITIES)

            self.firingPriorities(mode) = firingPriority;
        end

        function self = setFiringWeights(self, mode, firingWeight)
            % SELF = SETFIRINGWEIGHTS(MODE, FIRINGWEIGHTS)

            self.firingWeights(mode) = firingWeight;
        end

        function self = setFiringOutcome(self, mode, class, node, firingOutcome)
            % SELF = SETFIRINGOUTCOMES(MODE, NODE, CLASS,FIRINGOUTCOME)
            ind = self.model.getNodeIndex(node);
            self.firingOutcomes{mode}(ind,class) = firingOutcome;
        end

        function self = setDistribution(self, mode, distribution)
            self.distributions{mode} = distribution;
        end

        function self = setFiringRateDependence(self, mode, g)
            % SELF = SETFIRINGRATEDEPENDENCE(MODE, G)
            % Marking-dependent firing-rate multiplier for a timed mode.
            % G is a function handle g(m) of the input-place marking matrix m
            % (shaped like enablingConditions{mode}, nnodes x nclasses),
            % returning a positive scalar. The effective firing rate of an
            % enabled binding becomes rate_base(mode)*g(m). Empty G restores
            % the unit (marking-independent) multiplier.
            %
            % Exact only for memoryless firing: the mode must be TIMED (not
            % IMMEDIATE, which uses firingWeights) and exponentially
            % distributed. This mirrors the PS/FCFS-only restriction on
            % station load/class dependence.
            if ~isempty(g) && ~isa(g,'function_handle')
                line_error(mfilename,'Firing-rate dependence must be a function handle g(m) or empty.');
            end
            if ~isempty(g)
                if self.timingStrategies(mode) == TimingStrategy.IMMEDIATE
                    line_error(mfilename,'Firing-rate dependence is not supported for IMMEDIATE modes; use setFiringWeights for immediate transitions.');
                end
                if ~isa(self.distributions{mode},'Exp')
                    line_error(mfilename,'Firing-rate dependence is supported only for exponentially distributed (memoryless) firing; set an Exp distribution on the mode first.');
                end
            end
            self.firingRateDependence{mode} = g;
        end

        function nmodes = getNumberOfModes(self) 
            nmodes = length(self.modeNames);
        end

        function modes = getModes(self) 
            modes = self.modes;
        end

        function [map,mu,phi] = getServiceRates(self)
            % [PH,MU,PHI] = GETPHSERVICERATES()

            map = cell(1,self.getNumberOfModes);
            mu = cell(1,self.getNumberOfModes);
            phi = cell(1,self.getNumberOfModes);

            for r=1:self.getNumberOfModes
                switch class(self.distributions{r})
                    case {'Replayer','Trace'}
                        aph = self.distributions{r}.fitAPH;
                        map{r} = aph.getProcess();
                        mu{r} = aph.getMu;
                        phi{r} = aph.getPhi;
                    case {'Exp','Coxian','Erlang','HyperExp','Markovian','APH','MAP'}
                        map{r} = self.distributions{r}.getProcess();
                        mu{r} = self.distributions{r}.getMu;
                        phi{r} = self.distributions{r}.getPhi;
                    case 'MMPP2'
                        map{r} = self.distributions{r}.getProcess();
                        mu{r} = self.distributions{r}.getMu;
                        phi{r} = self.distributions{r}.getPhi;
                    case {'Det','Uniform','Pareto','Gamma','Weibull','Lognormal','Expolynomial'}
                        map{r} = self.distributions{r}.getProcess();
                        mu{r} = [self.distributions{r}.getRate];
                        phi{r} = [1];
                    otherwise
                        map{r}  = {[NaN],[NaN]};
                        mu{r}  = NaN;
                        phi{r}  = NaN;
                end
            end
        end
    end
end
