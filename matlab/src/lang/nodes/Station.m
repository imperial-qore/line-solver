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
        lcdScaling; % limited class-dependence scaling factors (product-form beta_{i,r})
        lcdScalingPeak; % peak (max) class-dependent rate scaling per class, used to normalize Util (T*S/peak)
        ljdScaling; % limited joint-dependence scaling factors (non-product-form eta_i)
        ljdScalingPeak; % peak (max) joint-dependent rate scaling per class, used to normalize Util (T*S/peak)
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
            self.patienceDistributions = [];
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedLoadDependence(self, alpha)
            % SETLIMITEDLOADDEPENDENCE(self, alpha)
            % alpha(ni) is the service rate scaling when there are ni>=1
            % jobs in the system
            self.lldScaling = alpha;
            self.invalidateStruct();
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedClassDependence(self, gamma, peakRatePerClass)
            % SETLIMITEDCLASSDEPENDENCE(self, gamma, peakRatePerClass)
            %
            % gamma(ni) is a function handle, where ni=[ni1,...,niR]
            % is the service rate scaling when there are nir jobs at
            % station i in class r.
            % peakRatePerClass is the peak (maximum) value of the rate
            % scaling per class (e.g. the effective number of servers). It is
            % REQUIRED and is used to normalize utilization as Util = T*S/peak,
            % matching the T*S/c convention of ordinary multiserver stations.
            % A scalar is broadcast to every class.
            if ~isa(gamma,'function_handle')
                line_error(mfilename, 'Class dependence must be specified through a function handle.');
            end
            if nargin < 3 || isempty(peakRatePerClass)
                line_error(mfilename, 'Class dependence requires an explicit peak rate: setClassDependence(beta, peakRatePerClass). Pass a scalar (identical peak for every class) or a per-class vector.');
            end
            if ~isnumeric(peakRatePerClass) || any(peakRatePerClass(:) <= 0)
                line_error(mfilename, 'peakRatePerClass must be a positive scalar or per-class vector.');
            end
            self.lcdScaling = gamma;
            self.lcdScalingPeak = peakRatePerClass(:)';
            self.invalidateStruct();
        end

        % don't expose to avoid accidental call without checking the queue
        % scheduling discipline
        function setLimitedJointDependence(self, eta, peakRatePerClass)
            % SETLIMITEDJOINTDEPENDENCE(self, eta, peakRatePerClass)
            %
            % eta(ni) is a function handle, where ni=[ni1,...,niR] is the
            % joint per-class population at the station. It returns either a
            % scalar service-rate scaling shared by every class, or a length-R
            % per-class vector. Unlike setLimitedClassDependence (beta_{i,r},
            % which must depend on the own-class marginal n_{i,r} and preserves
            % BCMP product form), eta may read the joint vector arbitrarily and
            % is therefore NON-product-form: solvers treat it as an
            % approximation with no exactness/uniqueness guarantee.
            % peakRatePerClass is REQUIRED and normalizes Util = T*S/peak; a
            % scalar is broadcast to every class.
            if ~isa(eta,'function_handle')
                line_error(mfilename, 'Joint dependence must be specified through a function handle.');
            end
            if nargin < 3 || isempty(peakRatePerClass)
                line_error(mfilename, 'Joint dependence requires an explicit peak rate: setJointDependence(eta, peakRatePerClass). Pass a scalar (identical peak for every class) or a per-class vector.');
            end
            if ~isnumeric(peakRatePerClass) || any(peakRatePerClass(:) <= 0)
                line_error(mfilename, 'peakRatePerClass must be a positive scalar or per-class vector.');
            end
            self.ljdScaling = eta;
            self.ljdScalingPeak = peakRatePerClass(:)';
            self.invalidateStruct();
        end

        function invalidateStruct(self)
            % INVALIDATESTRUCT(self)
            % Discard any cached NetworkStruct after a rate-scaling change.
            % Without this a setter called AFTER the first getStruct() is
            % silently dropped: the solver reads the stale sn and solves the
            % model unscaled, with no error.
            if ~isempty(self.model) && ismethod(self.model,'resetStruct')
                self.model.resetStruct();
            end
        end

    end

    methods

        function self = removeJobClass(self, jobclass)
            % SELF = REMOVEJOBCLASS(JOBCLASS)
            %
            % Drop the per-class capacity, drop rule and patience of JOBCLASS
            % on top of the routing configuration handled by Node.

            remaining = self.remainingClassIndexes(jobclass);
            removeJobClass@Node(self, jobclass);
            K = numel(remaining) + 1;
            if numel(self.classCap) == K
                self.classCap = self.classCap(remaining);
            end
            if numel(self.dropRule) == K
                self.dropRule = self.dropRule(remaining);
            end
            if numel(self.patienceDistributions) == K
                self.patienceDistributions = self.patienceDistributions(remaining);
            end
        end

        function self = setDropRule(self, class, drop)
            % SELF = SETDROPRULE(CLASS, DROPRULE)

            self.dropRule(class) = drop;
            self.invalidateStruct();
        end


        function setNumServers(self, value)
            % SETNUMSERVERS(VALUE)

            self.numberOfServers = value;
            self.invalidateStruct();
        end

        function setNumberOfServers(self, value)
            % SETNUMBEROFSERVERS(VALUE)

            self.numberOfServers = value;
            self.invalidateStruct();
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
            %
            % INVALIDATESTRUCT IS PART OF THE SETTER, not an optimization the
            % caller may skip: sn.cap and sn.classcap are DERIVED (refreshCapacity
            % folds this value together with classCap and the chain population),
            % so a cached struct does not see the new buffer. Without it, a
            % setCapacity called after the first getStruct() -- the ordinary
            % order when a model is built, inspected, then capped -- was silently
            % dropped and every sn-reading solver answered the UNBOUNDED model:
            % SolverCTMC returned the product-form 1.1475 jobs for a buffer of 1
            % while SolverFLD, whose gate reads the node objects instead, refused
            % the very same model as capacity-bound. see _kb/11-conventions-and-gotchas.md

            self.cap = value;
            self.invalidateStruct();
        end

        function setCap(self, value)
            % SETCAP(VALUE)
            % Alias for setCapacity() for backwards compatibility

            self.setCapacity(value);
        end

        function setClassCapacity(self, class, capacity)
            % SETCLASSCAPACITY(CLASS, CAPACITY)
            % Per-class buffer at this station, the station-level twin of
            % SETCHAINCAPACITY. Native Python (Station.set_class_capacity), C++
            % (network_builder set_class_capacity) and the JAR (setClassCap)
            % all carry it; MATLAB had it on Place only, so JSIM2LINE's import
            % of a per-class JSIMgraph capacity and every model that declares
            % one on a Queue died on an unrecognized method.
            %
            % SELF.CAP IS LEFT ALONE, unlike setChainCapacity, which sets EVERY
            % class in one call and can therefore total them. Setting one class
            % says nothing about the others, and REFRESHCAPACITY already reads
            % classCap(r) beside cap and takes the tighter of the two.

            % Resolve a JobClass object to its column index (mirrors Place.m);
            % indexing classCap by the object itself silently fails to store the
            % capacity, so the station stays unbounded and this is a no-op.
            if isa(class, 'JobClass')
                r = class.index;
            else
                r = class;
            end
            if ~(capacity > 0)
                line_error(mfilename, sprintf(['Class capacity must be positive (Inf for unbounded), got %g. ' ...
                    'A zero capacity is how refreshCapacity encodes a class the station does not serve, ' ...
                    'so it cannot also mean a buffer of size zero.'], capacity));
            end
            self.classCap(r) = capacity;
            self.invalidateStruct();
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
            % This one MATERIALIZES the struct itself (getStruct above) before
            % mutating classCap/cap, so the cache is guaranteed stale on return.
            self.invalidateStruct();
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
                        case {'Exp','Coxian','Erlang','HyperExp','Markovian','APH','PH','ME','CME','RAP'}
                            map{r} = self.input.sourceClasses{r}{end}.getProcess;
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case {'Det','Uniform','Pareto','Gamma','Lognormal','Weibull'}
                            map{r}  = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r}  = [1/self.input.sourceClasses{r}{end}.getMean];
                            phi{r}  = [1];
                        case {'Bernoulli','Binomial','Poisson'}
                            % see _kb/04-networkstruct.md (node/process construction notes) for rationale
                            map{r}  = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r}  = [self.input.sourceClasses{r}{end}.getRate];
                            phi{r}  = [1];
                        case 'Geometric'
                            % see _kb/04-networkstruct.md (node/process construction notes) for rationale
                            map{r}  = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r}  = [self.input.sourceClasses{r}{end}.getRate];
                            phi{r}  = [1];
                        case 'DMAP'
                            % Discrete-time (D0,D1): the pair already has MAP
                            % shape, so it reaches sn.proc verbatim. Without
                            % this arm the switch fell through and a DMAP
                            % source arrived with an EMPTY representation.
                            % mu/phi carry the per-slot event rate, which is
                            % all a non-slotted consumer can read from it.
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r} = 1 / self.input.sourceClasses{r}{end}.getMean;
                            phi{r} = 1;
                        case 'MMPP2'
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case 'MAP'
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r} = self.input.sourceClasses{r}{end}.getMu;
                            phi{r} = self.input.sourceClasses{r}{end}.getPhi;
                        case 'NHPP'
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mu{r} = self.input.sourceClasses{r}{end}.getTimeAverageRate();
                            phi{r} = 1;
                        case {'MAPt','PHt'}
                            % The schedule is the parameterisation; mu/phi carry
                            % the time-averaged nominal so the phase count is
                            % preserved, unlike NHPP which collapses to one phase.
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            if strcmp(class(self.input.sourceClasses{r}{end}), 'MAPt')
                                [D0bar, D1bar] = self.input.sourceClasses{r}{end}.getTimeAverageProcess();
                            else
                                [D0bar, D1bar] = self.input.sourceClasses{r}{end}.getTimeAverageProcessMAP();
                            end
                            mu{r} = -diag(D0bar);
                            phi{r} = sum(D1bar, 2) ./ mu{r};
                        case 'BMAP'
                            % {D0, D1, D_batch1, ..., D_batchK}: same layout
                            % as the JAR MatrixCell for batch arrivals
                            map{r} = self.input.sourceClasses{r}{end}.getProcess();
                            mapAggr = self.input.sourceClasses{r}{end}.toMAP;
                            mu{r} = mapAggr.getMu;
                            phi{r} = mapAggr.getPhi;
                        case 'MarkedMAP'
                            % see _kb/04-networkstruct.md (node/process construction notes) for rationale
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
            
            % SIZED BY THE MODEL'S CLASS COUNT, NOT BY THE CELL. Sizing it from
            % the cell returns SHORT output whenever a station's per-class slots
            % were never padded to nclasses, and every caller then walks r=1:K
            % over it and throws MATLAB's own "Index exceeds the number of array
            % elements" -- naming neither the station nor the class, and burying
            % the real fault (a class with no service configured there). The
            % absent slots fall into the isempty arm below, which is what the
            % padding would have produced anyway.
            nclasses = size(self.server.serviceProcess,2);
            if ~isempty(self.model) && ismethod(self.model,'getNumberOfClasses')
                nclasses = max(nclasses, self.model.getNumberOfClasses());
            end
            map = cell(1,nclasses);
            mu = cell(1,nclasses);
            phi = cell(1,nclasses);
            for r=1:nclasses
                serviceProcess_r = [];
                if r <= numel(self.server.serviceProcess)
                    serviceProcess_r = self.server.serviceProcess{r};
                end
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
                        case {'Bernoulli','Binomial','Poisson'}
                            % Counting distributions used as a service time; see
                            % the matching arm in getSourceRates above.
                            map{r}  = serviceProcess_r{end}.getProcess();
                            mu{r}  = [serviceProcess_r{end}.getRate];
                            phi{r}  = [1];
                        case 'Geometric'
                            % Lattice-valued service time supported on
                            % {1,2,...}; see the matching arm in
                            % getSourceRates above.
                            map{r}  = serviceProcess_r{end}.getProcess();
                            mu{r}  = [serviceProcess_r{end}.getRate];
                            phi{r}  = [1];
                        case {'Replayer', 'Trace'}
                            aph = serviceProcess_r{end}.fitAPH;
                            map{r} = aph.getProcess();
                            mu{r} = aph.getMu;
                            phi{r} = aph.getPhi;
                        case {'Exp','Coxian','Erlang','HyperExp','Markovian','APH','MAP','PH','ME','CME','RAP'}
                            map{r} = serviceProcess_r{end}.getProcess();
                            mu{r} = serviceProcess_r{end}.getMu;
                            phi{r} = serviceProcess_r{end}.getPhi;
                        case 'NHPP'
                            % Mirrors the Source arm above. Without this case the
                            % switch falls through with map/mu/phi unassigned, so a
                            % time-varying service rate reached sn.proc empty.
                            map{r} = serviceProcess_r{end}.getProcess();
                            mu{r} = serviceProcess_r{end}.getTimeAverageRate();
                            phi{r} = 1;
                        case {'MAPt','PHt'}
                            % Mirrors the Source arm: the nominal preserves the
                            % phase count that the schedule modulates.
                            map{r} = serviceProcess_r{end}.getProcess();
                            if strcmp(class(serviceProcess_r{end}), 'MAPt')
                                [D0bar, D1bar] = serviceProcess_r{end}.getTimeAverageProcess();
                            else
                                [D0bar, D1bar] = serviceProcess_r{end}.getTimeAverageProcessMAP();
                            end
                            mu{r} = -diag(D0bar);
                            phi{r} = sum(D1bar, 2) ./ mu{r};
                        case 'MMPP2'
                            map{r} = serviceProcess_r{end}.getProcess();
                            mu{r} = serviceProcess_r{end}.getMu;
                            phi{r} = serviceProcess_r{end}.getPhi;
                        case 'DMAP'
                            % Mirrors the Source arm: the discrete (D0,D1)
                            % reaches sn.proc verbatim, mu/phi hold the
                            % per-slot event rate.
                            map{r} = serviceProcess_r{end}.getProcess();
                            mu{r} = 1 / serviceProcess_r{end}.getMean;
                            phi{r} = 1;
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
