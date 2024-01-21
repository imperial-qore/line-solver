classdef Transition < Node
    % A service station with queueing
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        enablingConditions;
        inhibitingConditions;
        modeNames;
        numbersOfServers;
        timingStrategies;
        distributions;
        firingPriorities;
        firingWeights;
        firingOutcomes;
        % schedStrategy;
        % schedStrategyPar;
        nmodes;
        cap;
    end
    
    methods
        function self = Transition(model,name)
            % TRANSITION(MODEL, NAME)
            
            self@Node(name);
            classes = model.classes;
            self.input = Enabling(classes);
            self.output = Firing(classes);
            self.cap = Inf; % Compatible with other nodes
            
            self.setModel(model);
            self.model.addNode(self);

            self.server = Timing();

            self.enablingConditions = [];
            self.inhibitingConditions = [];
            self.modeNames = {};
            self.numbersOfServers = [];
            self.timingStrategies = [];
            self.distributions = [];
            self.firingPriorities = [];
            self.firingWeights = [];
            self.firingOutcomes = [];
            self.nmodes = 0;
        end
   
        function self = init(self)
            % SELF = INIT()
            
            nclasses = length(self.model.classes);
            nnodes = length(self.model.nodes);
            self.nmodes = length(self.modeNames);

            self.enablingConditions = cell(self.nmodes);
            self.inhibitingConditions = cell(self.nmodes);
            self.firingOutcomes = cell(self.nmodes);
            for m=1:self.nmodes
                self.enablingConditions{m} = zeros(nnodes,nclasses);            
                self.inhibitingConditions{m} = zeros(nnodes,nclasses);            
                self.firingOutcomes{m} = zeros(nnodes,nclasses);                            
            end
            self.numbersOfServers = Inf(1,self.nmodes);
            self.timingStrategies = repmat(TimingStrategy.ID_TIMED,1,self.nmodes);
            self.firingWeights = ones(1,self.nmodes);
            self.firingPriorities = ones(1,self.nmodes);
            self.distributions = cell(1, self.nmodes);
            self.distributions(:) = {Exp(1)};
        end
        
        function self = addMode(self, modeName)
            self.modeNames{end+1} = modeName;
            self.enablingConditions(:,end+1) = 0;
            self.inhibitingConditions(:,end+1) = Inf;
            self.numbersOfServers(end+1) = Inf;
            self.timingStrategies(end+1) = TimingStrategy.ID_TIMED;
            self.firingWeights(end+1) = 1;
            self.firingPriorities(end+1) = 1;
            self.distributions{end+1} = Exp(1);
            self.firingOutcomes(:,end+1) = 0;
            self.nmodes = self.nmodes + 1;            
        end

        function self = setEnablingConditions(self, mode, class, node, enablingCondition)
            % SELF = SETENABLINGCONDITIONS(MODE, CLASS, NODE, ENABLINGCONDITIONS)
            
            nnodes = length(self.model.nodes);
            if isa(node, 'Place')
                node = self.model.getNodeIndex(node.name);
                self.enablingConditions{mode}(node,class) = enablingCondition;
            elseif isnumeric(node) && node <= nnodes && isa(self.model.nodes{node}, 'Place')
                self.enablingConditions{mode}(node,class) = enablingCondition;
            else
                error('Node must be a Place node or index of a Place node.');
            end
        end

        function self = setInhibitingConditions(self, mode, class, node, inhibitingCondition)
            % SELF = SETINHIBITINGCONDITIONS(MODE, CLASS, NODE, INHIBITINGCONDITIONS)
            
            nnodes = length(self.model.nodes);
            if isa(node, 'Place')
                node = self.model.getNodeIndex(node.name);
                self.inhibitingConditions{mode}(node,class) = inhibitingCondition;
            elseif isa(node, 'double') && node <= nnodes && isa(self.model.nodes{node}, 'Place')
                self.inhibitingConditions{mode}(node,class) = inhibitingCondition;
            else
                error('Node must be a Place node or index of a Place node.');
            end
        end

        function self = setModeNames(self, mode, modeName)
            % SELF = SETMODENAMES(MODE, MODENAMES)
            
            self.modeNames{mode} = modeName;
        end

        function self = setNumberOfServers(self, mode, numberOfServers)
            % SELF = SETNUMBEROFSERVERS(MODE, NUMOFSERVERS)
            
            self.numbersOfServers(mode) = numberOfServers;
        end

        function self = setTimingStrategy(self, mode, timingStrategy)
            % SELF = SETTIMINGSTRATEGY(MODE, TIMINGSTRATEGY)

            self.timingStrategies(mode) = TimingStrategy.toId(timingStrategy);
        end
        
        function self = setFiringPriorities(self, mode, firingPriority)
            % SELF = SETFIRINGPRIORITIES(MODE, FIRINGPRIORITIES)
            
            self.firingPriorities(mode) = firingPriority;
        end

        function self = setFiringWeights(self, mode, firingWeight)
            % SELF = SETFIRINGWEIGHTS(MODE, FIRINGWEIGHTS)
            
            self.firingWeights(mode) = firingWeight;
        end

        function self = setFiringOutcome(self, class, mode, node, firingOutcome)
            % SELF = SETFIRINGOUTCOMES(MODE, CLASS, NODE, FIRINGOUTCOME)

            nnodes = length(self.model.nodes);
            if isa(node, 'Node')
                node = self.model.getNodeIndex(node.name);
                self.firingOutcomes{mode}(node,class) = firingOutcome;
            elseif isa(node, 'double') && node <= nnodes && isa(self.model.nodes{node}, 'Node')
                self.firingOutcomes{mode}(node,class) = firingOutcome;
            else
                error('Node is not valid.');
            end
        end

        function self = setDistribution(self, mode, distribution)
            self.distributions{mode} = distribution;
        end

        function bool = isTransition(self)
            % BOOL = ISTRANSITION()
            
            bool = isa(self,'Transition');
        end
        
        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()
            
            sections = {self.input, self.server, self.output};
        end


        function [map,mu,phi] = getServiceRates(self)
            % [PH,MU,PHI] = GETPHSERVICERATES()
            
            nmodes = self.nmodes;
            map = cell(1,nmodes);
            mu = cell(1,nmodes);
            phi = cell(1,nmodes);

            for r=1:nmodes
                switch class(self.distributions{r})
                    case 'Replayer'
                        aph = self.distributions{r}.fitAPH;
                        map{r} = aph.getRepres();
                        mu{r} = aph.getMu;
                        phi{r} = aph.getPhi;
                    case {'Exp','Coxian','Erlang','HyperExp','MarkovianDistribution','APH','MAP'}
                        map{r} = self.distributions{r}.getRepres();
                        mu{r} = self.distributions{r}.getMu;
                        phi{r} = self.distributions{r}.getPhi;
                    case 'MMPP2'
                        map{r} = self.distributions{r}.getRepres();
                        mu{r} = self.distributions{r}.getMu;
                        phi{r} = self.distributions{r}.getPhi;
                    otherwise
                        map{r}  = {[NaN],[NaN]};
                        mu{r}  = NaN;
                        phi{r}  = NaN;
                end
            end
        end
    end
end
