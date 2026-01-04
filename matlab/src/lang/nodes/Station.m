classdef Station < StatefulNode
    % An abstract class for nodes where jobs station
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        numberOfServers;
        cap;
        dropRule;
        classCap;
        lldScaling; % limited load-dependence scaling factors
        lcdScaling; % limited class-dependence scaling factors
        ljdScaling; % limited joint-dependence scaling factors (linearized)
        ljdCutoffs; % per-class cutoffs for joint dependence [N1, N2, ..., NR]
        ljcdScaling; % limited joint-class-dependence scaling factors (cell array, per class)
        ljcdCutoffs; % per-class cutoffs for joint class dependence [N1, N2, ..., NR]
        stationIndex;
        patienceDistributions; % per-class patience distributions (cell array indexed by class)
    end

    methods(Hidden)
        %Constructor
        function self = Station(name)
            % SELF = STATION(NAME)

            self@StatefulNode(name);
            self.cap = Inf;
            self.classCap = [];
            self.lldScaling = [];
            self.ljdScaling = [];
            self.ljdCutoffs = [];
            self.ljcdScaling = {};
            self.ljcdCutoffs = [];
            self.patienceDistributions = [];
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedLoadDependence(self, alpha)
            % SETLIMITEDLOADDEPENDENCE(self, alpha)
            % alpha(ni) is the service rate scaling when there are ni>=1
            % jobs in the system
            self.lldScaling = alpha;
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedClassDependence(self, gamma)
            % SETLIMITEDCLASSDEPENDENCE(self, gamma)
            %
            % gamma(ni) is a function handle, where ni=[ni1,...,niR]
            % is the service rate scaling when there are nir jobs at
            % station i in class r
            if isa(gamma,'function_handle')
                self.lcdScaling = gamma;
            else
                line_error(mfilename, 'Class dependence must be specified through a function handle.');
            end
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedJointDependence(self, scalingTable, cutoffs)
            % SETLIMITEDJOINTDEPENDENCE(self, scalingTable, cutoffs)
            %
            % scalingTable is a linearized vector of scaling factors
            % indexed by: idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...
            % cutoffs is a vector [N1, N2, ..., NR] of per-class cutoffs
            self.ljdScaling = scalingTable(:)';
            self.ljdCutoffs = cutoffs(:)';
        end

        function [scalingTable, cutoffs] = getLimitedJointDependence(self)
            % [SCALINGTABLE, CUTOFFS] = GETLIMITEDJOINTDEPENDENCE(self)
            %
            % Returns the joint-dependence scaling table and cutoffs
            scalingTable = self.ljdScaling;
            cutoffs = self.ljdCutoffs;
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedJointClassDependence(self, scalingTables, cutoffs)
            % SETLIMITEDJOINTCLASSDEPENDENCE(self, scalingTables, cutoffs)
            %
            % scalingTables is a cell array {1,K} where each cell contains
            % a linearized vector of scaling factors for that class
            % indexed by: idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...
            % cutoffs is a vector [N1, N2, ..., NR] of per-class cutoffs
            self.ljcdScaling = scalingTables;
            self.ljcdCutoffs = cutoffs(:)';
        end

        function [scalingTables, cutoffs] = getLimitedJointClassDependence(self)
            % [SCALINGTABLES, CUTOFFS] = GETLIMITEDJOINTCLASSDEPENDENCE(self)
            %
            % Returns the joint-class-dependence scaling tables and cutoffs
            scalingTables = self.ljcdScaling;
            cutoffs = self.ljcdCutoffs;
        end

    end

    methods

        function self = setDropRule(self, class, drop)
            % SELF = SETDROPRULE(CLASS, DROPRULE)

            self.dropRule(class) = drop;
        end


        function setNumServers(self, value)
            % SETNUMSERVERS(VALUE)

            self.numberOfServers = value;
        end

        function setNumberOfServers(self, value)
            % SETNUMBEROFSERVERS(VALUE)

            self.numberOfServers = value;
        end

        function value = getNumServers(self)
            % VALUE = GETNUMSERVERS()

            value = self.numberOfServers;
        end

        function value = getNumberOfServers(self)
            % VALUE = GETNUMBEROFSERVERS()

            value = self.numberOfServers;
        end

        function setCapacity(self, value)
            % SETCAPACITY(VALUE)

            self.cap = value;
        end

        function setCap(self, value)
            % SETCAP(VALUE)
            % Alias for setCapacity() for backwards compatibility

            self.setCapacity(value);
        end

        function setChainCapacity(self, values)
            % SETCHAINCAPACITY(VALUES)

            sn = self.model.getStruct;
            if numel(values) ~= sn.nchains
                line_error(mfilename,'The method requires in input a capacity value for each chain.');
            end
            for c = 1:sn.nchains
                inchain = sn.inchain{c};
                for r = inchain
                    if ~self.isServiceDisabled(r)
                        self.classCap(r) = values(c);
                    else
                        self.classCap(r) = Inf;
                    end
                end
            end
            self.cap = min(sum(self.classCap(self.classCap>0)), self.cap);
        end


        function isD = isServiceDefined(self, class)
            K = size(self.model.getClasses(),2);
            isD = true(1, K);
            switch self.server.className
                case 'ServiceTunnel'
                    %noop
                otherwise
                    for r=1:K
                        if isempty(self.server.serviceProcess{1,r})
                            isD(r) = false;
                        end
                    end
            end
        end

        function isD = isServiceDisabled(self, class)
            % ISD = ISSERVICEDISABLED(CLASS)
            if nargin>=2
                switch self.server.className
                    case 'ServiceTunnel'
                        isD = false;
                    otherwise
                        isD = self.server.serviceProcess{1,class}{end}.isDisabled();
                end
            else
                K = size(self.model.getClasses(),2);
                isD = false(1, K);
                switch self.server.className
                    case 'ServiceTunnel'
                        %noop
                    otherwise
                        %isD = cellfun(@(sp) sp{end}.isDisabled, self.server.serviceProcess);
                        for r=1:K
                            isD(r) = self.server.serviceProcess{1,r}{end}.isDisabled();
                        end
                end
            end
        end

        function isI = isServiceImmediate(self, class)
            % ISI = ISSERVICEIMMEDIATE(CLASS)

            isI = self.server.serviceProcess{1,class}{end}.isImmediate();
        end

        function R = getNumberOfServiceClasses(self)
            % R = GETNUMBEROFSERVICECLASSES()

            R = size(self.server.serviceProcess,2);
        end

        function [p] = getSelfLoopProbabilities(self)
            % [P] = GETSELFLOOPPROBABILITIES()

            R = getNumberOfServiceClasses(self);
            p = zeros(1,R);
            for k=1:R
                nOutLinks = length(self.output.outputStrategy{k}{end});
                switch RoutingStrategy.toText(self.output.outputStrategy{k}{2})
                    case 'Random'
                        p(k) = 1 / nOutLinks;
                    case RoutingStrategy.PROB
                        for t=1:nOutLinks % for all outgoing links
                            if strcmp(self.output.outputStrategy{k}{end}{t}{1}.name, self.name)
                                p(k) = self.output.outputStrategy{k}{end}{t}{2};
                                break
                            end
                        end
                end
            end
        end

        function [map, mu, phi] = getSourceRates(self)
            % [PH,MU,PHI] = GETSOURCERATES()

            nclasses = size(self.input.sourceClasses,2);
            map = cell(1,nclasses);
            mu = cell(1,nclasses);
            phi = cell(1,nclasses);
            for r=1:nclasses
                if isempty(self.input.sourceClasses{r})
                    self.input.sourceClasses{r} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                elseif ~self.input.sourceClasses{r}{end}.isDisabled()
                    switch class(self.input.sourceClasses{r}{end})
                        case {'Replayer', 'Trace'}
                            aph = self.input.sourceClasses{r}{end}.fitAPH;
                            map{r} = aph.getProcess();
                            mu{r} = aph.getMu;
                            phi{r} = aph.getPhi;
                        case {'Exp','Coxian','Erlang','HyperExp','Markovian','APH','PH'}
                            map{r} = self.input.sourceClasses{r}{end}.getProcess;
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case {'Det','Uniform','Pareto','Gamma','Lognormal','Weibull'}
                            map{r}  = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r}  = [1/self.input.sourceClasses{r}{end}.getMean];
                            phi{r}  = [1];
                        case 'MMPP2'
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case 'MAP'
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        otherwise
                            % leave everything empty
                    end
                else
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                end
            end
        end

        function [map,mu,phi] = getServiceRates(self)
            % [PH,MU,PHI] = GETSERVICERATES()
            
            nclasses = size(self.server.serviceProcess,2);
            map = cell(1,nclasses);
            mu = cell(1,nclasses);
            phi = cell(1,nclasses);
            for r=1:nclasses
                serviceProcess_r = self.server.serviceProcess{r};
                if isempty(serviceProcess_r)
                    serviceProcess_r = {[],ServiceStrategy.LI,Disabled.getInstance()};
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                elseif serviceProcess_r{end}.isImmediate()
                    map{r}  = {[-GlobalConstants.Immediate],[GlobalConstants.Immediate]};
                    mu{r}  = [GlobalConstants.Immediate];
                    phi{r}  = [1];
                elseif ~serviceProcess_r{end}.isDisabled()
                    switch class(serviceProcess_r{end})
                        case {'Det','Uniform','Pareto','Gamma','Weibull','Lognormal'}
                            map{r}  = serviceProcess_r{end}.getProcess();
                            mu{r}  = [serviceProcess_r{end}.getRate];
                            phi{r}  = [1];
                        case {'Replayer', 'Trace'}
                            aph = serviceProcess_r{end}.fitAPH;
                            map{r} = aph.getProcess();
                            mu{r} = aph.getMu;
                            phi{r} = aph.getPhi;
                        case {'Exp','Coxian','Erlang','HyperExp','Markovian','APH','MAP','PH'}
                            map{r} = serviceProcess_r{end}.getProcess();
                            mu{r} = serviceProcess_r{end}.getMu;
                            phi{r} = serviceProcess_r{end}.getPhi;
                        case 'MMPP2'
                            map{r} = serviceProcess_r{end}.getProcess();
                            mu{r} = serviceProcess_r{end}.getMu;
                            phi{r} = serviceProcess_r{end}.getPhi;
                    end
                else
                    map{r}  = {[NaN],[NaN]};
                    mu{r}  = NaN;
                    phi{r}  = NaN;
                end
            end
        end

        function summary(self)
            % SUMMARY()

            line_printf('\nNode: <strong>%s</strong>',self.getName);
            line_printf('\nScheduling: %s',self.schedStrategy);
            line_printf('\nNumber of Servers: %d',self.numberOfServers);
            for r=1:length(self.output.outputStrategy)
                classes = self.model.getClasses();
                line_printf('\nRouting %s: %s',classes{r}.name,self.output.outputStrategy{r}{2});
            end
            %            self.input.summary;
            %            self.server.summary;
            %            self.output.summary;
        end

    end
end
