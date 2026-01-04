classdef MNetwork < Model 
    % An extended queueing network model.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Access=private)
        enableChecks;
        hasState;
        logPath;
        usedFeatures; % structure of booleans listing the used classes
        % it must be accessed via getUsedLangFeatures that updates
        % the Distribution classes dynamically
    end

    properties (Hidden)
        obj; % empty
        sn;
        csMatrix;
        hasStruct;
        allowReplace;
    end

    properties (Hidden)
        handles;
        sourceidx; % cached value
        sinkidx; % cached value
    end

    properties
        classes;
        items;
        stations;
        nodes;
        connections;
        regions;
    end

    methods % get methods
        nodes = getNodes(self)
        sa = getStruct(self, structType, wantState) % get abritrary representation
        used = getUsedLangFeatures(self) % get used features
        ft = getForkJoins(self, rt) % get fork-join pairs
        [chainsObj,chainsMatrix] = getChains(self, rt) % get chain table
        [rt,rtNodes,connections,chains,rtNodesByClass,rtNodesByStation] = getRoutingMatrix(self, arvRates) % get routing matrix
        [Q,U,R,T,A,W] = getAvgHandles(self)
        Q = getAvgQLenHandles(self);
        U = getAvgUtilHandles(self);
        R = getAvgRespTHandles(self);
        T = getAvgTputHandles(self);
        A = getAvgArvRHandles(self);
        W = getAvgResidTHandles(self);
        [Qt,Ut,Tt] = getTranHandles(self)
        connections = getConnectionMatrix(self);
        [absorbingStations, absorbingIdxs] = getAbsorbingStations(self) % get absorbing stations
        info = getReducibilityInfo(self) % get reducibility analysis
    end

    methods % ergodicity analysis methods
        [isErg, info] = isRoutingErgodic(self, P) % check if routing is ergodic
        P = makeErgodic(self, targetNode) % generate ergodic routing matrix
    end

    methods % link, reset, refresh methods
        self = link(self, P)
        self = relink(self, P)
        [loggerBefore,loggerAfter] = linkAndLog(self, nodes, classes, P, wantLogger, logPath)
        sanitize(self);

        reset(self, resetState)
        resetHandles(self)
        resetModel(self, resetState)
        nodes = resetNetwork(self, deleteCSnodes)
        resetStruct(self)

        refreshStruct(self, hard);
        [rates, scv, hasRateChanged, hasSCVChanged] = refreshRates(self, statSet, classSet);
        [ph, mu, phi, phases] = refresProcessPhases(self, statSet, classSet);
        proctypes = refreshProcessTypes(self);
        [rt, rtfun, rtnodes] = refreshRoutingMatrix(self, rates);
        [lt] = refreshLST(self, statSet, classSet);
        sync = refreshSync(self);
        classprio = refreshPriorities(self);
        [sched, schedparam] = refreshScheduling(self);
        [rates, scv, mu, phi, phases] = refreshProcesses(self, statSet, classSet);
        [chains, visits, rt] = refreshChains(self, propagate)
        [cap, classcap] = refreshCapacity(self);
        nvars = refreshLocalVars(self);
    end

    % PUBLIC METHODS
    methods (Access=public)

        %Constructor
        function self = MNetwork(modelName, varargin)
            % SELF = NETWORK(MODELNAME)
            self@Model(modelName);
            self.nodes = {};
            self.stations = {};
            self.classes = {};
            self.connections = [];
            initUsedFeatures(self);
            self.sn = [];
            self.hasState = false;
            self.logPath = '';
            self.items = {};
            self.regions = {};
            self.sourceidx = [];
            self.sinkidx = [];
            self.setChecks(true);
            self.hasStruct = false;
            self.allowReplace = false;
        end

        setInitialized(self, bool);

        function self = setChecks(self, bool)
            self.enableChecks = bool;
        end

        P = getLinkedRoutingMatrix(self)

        function logPath = getLogPath(self)
            % LOGPATH = GETLOGPATH()

            logPath = self.logPath;
        end

        function setLogPath(self, logPath)
            % SETLOGPATH(LOGPATH)

            self.logPath = logPath;
        end

        bool = hasInitState(self)

        function [M,R] = getSize(self)
            % [M,R] = GETSIZE()

            M = self.getNumberOfNodes;
            R = self.getNumberOfClasses;
        end

        function bool = hasOpenClasses(self)
            % BOOL = HASOPENCLASSES()

            bool = any(isinf(getNumberOfJobs(self)));
        end

        function bool = hasClassSwitching(self)
            % BOOL = HASCLASSSWITCHING()

            bool = any(cellfun(@(c) isa(c,'ClassSwitch'), self.nodes));
        end

        function bool = hasFork(self)
            % BOOL = HASFORK()

            bool = any(cellfun(@(c) isa(c,'Fork'), self.nodes));
        end

        function bool = hasJoin(self)
            % BOOL = HASJOIN()

            bool = any(cellfun(@(c) isa(c,'Join'), self.nodes));
        end

        function bool = hasClosedClasses(self)
            % BOOL = HASCLOSEDCLASSES()

            bool = any(isfinite(getNumberOfJobs(self)));
        end

        function bool = isMatlabNative(self)
            % BOOL = ISMATLABNATIVE()
            %
            % Returns true for MATLAB native (MNetwork) implementation
            
            bool = true;
        end

        function bool = isJavaNative(self)
            % BOOL = ISJAVANATIVE()
            %
            % Returns false for MATLAB native (MNetwork) implementation
            
            bool = false;
        end

        function index = getIndexOpenClasses(self)
            % INDEX = GETINDEXOPENCLASSES()

            index = find(isinf(getNumberOfJobs(self)))';
        end

        function index = getIndexClosedClasses(self)
            % INDEX = GETINDEXCLOSEDCLASSES()

            index = find(isfinite(getNumberOfJobs(self)))';
        end

        function classes = getClasses(self)
            classes= self.classes;
        end

        chain = getClassChain(self, className)
        c = getClassChainIndex(self, className)
        classNames = getClassNames(self)

        nodeNames = getNodeNames(self)
        nodeTypes = getNodeTypes(self)

        P = initRoutingMatrix(self)

        ind = getNodeIndex(self, name)
        lldScaling = getLimitedLoadDependence(self)
        lcdScaling = getLimitedClassDependence(self)

        function stationIndex = getStationIndex(self, name)
            % STATIONINDEX = GETSTATIONINDEX(NAME)

            if isa(name,'Node')
                node = name;
                name = node.getName();
            end
            stationIndex = find(cellfun(@(c) strcmp(c,name),self.getStationNames));
        end

        function statefulIndex = getStatefulNodeIndex(self, name)
            % STATEFULINDEX = GETSTATEFULNODEINDEX(NAME)

            if isa(name,'Node')
                node = name;
                name = node.getName();
            end
            statefulIndex = find(cellfun(@(c) strcmp(c,name),self.getStatefulNodeNames));
        end

        function classIndex = getClassIndex(self, name)
            % CLASSINDEX = GETCLASSINDEX(NAME)
            if isa(name,'JobClass')
                jobclass = name;
                name = jobclass.getName();
            end
            classIndex = find(cellfun(@(c) strcmp(c,name),self.getClassNames));
        end

        function stationnames = getStationNames(self)
            % STATIONNAMES = GETSTATIONNAMES()

            if self.hasStruct
                nodenames = self.sn.nodenames;
                isstation = self.sn.isstation;
                stationnames = {nodenames{isstation}}';
            else
                stationnames = {};
                for i=self.getStationIndexes
                    stationnames{end+1,1} = self.nodes{i}.name;
                end
            end
        end

        function nodes = getNodeByName(self, name)
            % NODES = GETNODEBYNAME(SELF, NAME)
            idx = findstring(self.getNodeNames,name);
            if idx > 0
                nodes = self.nodes{idx};
            else
                nodes = NaN;
            end
        end

        function station = getStationByName(self, name)
            % STATION = GETSTATIONBYNAME(SELF, NAME)
            idx = findstring(self.getStationNames,name);
            if idx > 0
                station = self.stations{idx};
            else
                station = NaN;
            end
        end

        function class = getClassByName(self, name)
            % CLASS = GETCLASSBYNAME(SELF, NAME)
            idx = findstring(self.getClassNames,name);
            if idx > 0
                class = self.classes{idx};
            else
                class = NaN;
            end
        end

        function nodes = getNodeByIndex(self, idx)
            % NODES = GETNODEBYINDEX(SELF, NAME)
            if idx > 0
                nodes = self.nodes{idx};
            else
                nodes = NaN;
            end
        end

        function station = getStationByIndex(self, idx)
            % STATION = GETSTATIONBYINDEX(SELF, NAME)
            if idx > 0
                station = self.stations{idx};
            else
                station = NaN;
            end
        end

        function class = getClassByIndex(self, idx)
            % CLASS = GETCLASSBYINDEX(SELF, NAME)
            if idx > 0
                class = self.classes{idx};
            else
                class = NaN;
            end
        end

        function [infGen, eventFilt, ev] =  getGenerator(self, varargin)
            line_warning(mfilename,'Results will not be cached. Use SolverCTMC(model,...).getGenerator(...) instead.\n');
            [infGen, eventFilt, ev] = SolverCTMC(self).getGenerator(varargin{:});
        end

        function [stateSpace,nodeStateSpace] = getStateSpace(self, varargin)
            line_warning(mfilename,'Results will not be cached. Use SolverCTMC(model,...).getStateSpace(...) instead.\n');
            [stateSpace,nodeStateSpace] = SolverCTMC(self).getStateSpace(varargin{:});
        end

        function summary(self)
            % SUMMARY()
            for i=1:self.getNumberOfClasses
                self.classes{i}.summary();
                if i<self.getNumberOfClasses
                    line_printf('\n');
                end
            end
            for i=1:self.getNumberOfNodes
                self.nodes{i}.summary();
            end
            line_printf('\n<strong>Routing matrix</strong>:');
            self.printRoutingMatrix
            line_printf('\n<strong>Product-form parameters</strong>:');
            [arvRates,servDemands,nJobs,thinkTimes,ldScalings,nServers]= getProductFormParameters(self);
            line_printf('\nArrival rates: %s',mat2str(arvRates,6));
            line_printf('\nService demands: %s',mat2str(servDemands,6));
            line_printf('\nNumber of jobs: %s',mat2str(nJobs));
            line_printf('\nThink times: %s',mat2str(thinkTimes,6));
            line_printf('\nLoad-dependent scalings: %s',mat2str(ldScalings,6));
            line_printf('\nNumber of servers: %s\n',mat2str(nServers));
            line_printf('\n<strong>Chains</strong>:');
        end

        function [D,Z] = getDemands(self)
            % [D,Z]= GETDEMANDS()
            %
            % Outputs:
            % D: service demands
            % Z: think times

            [~,D,~,Z,~,~] = sn_get_product_form_params(self.getStruct);
        end

        function [lambda,D,N,Z,mu,S]= getProductFormParameters(self)
            % [LAMBDA,D,N,Z,MU,S]= GETPRODUCTFORMPARAMETERS()
            %
            % Outputs:
            % LAMBDA: arrival rates for open classes
            % D: service demands
            % N: population vector for closed classes
            % Z: think times
            % mu: load-dependent rates
            % S: number of servers

            % mu also returns max(S) elements after population |N| as this is
            % required by MVALDMX

            [lambda,D,N,Z,mu,S] = sn_get_product_form_params(self.getStruct);
        end

        function [lambda,D,N,Z,mu,S]= getProductFormChainParameters(self)
            % [LAMBDA,D,N,Z,MU,S]= GETPRODUCTFORMCHAINPARAMETERS()
            %
            % Outputs:
            % LAMBDA: arrival rates for open classes
            % D: service demands
            % N: population vector for closed classes
            % Z: think times
            % mu: load-dependent rates
            % S: number of servers

            % mu also returns max(S) elements after population |N| as this is
            % required by MVALDMX

            [lambda,D,N,Z,mu,S]=  sn_get_product_form_chain_params(self.getStruct);
        end

        function statefulnodes = getStatefulNodes(self)
            statefulnodes = {};
            for i=1:self.getNumberOfNodes
                if self.nodes{i}.isStateful
                    statefulnodes{end+1,1} = self.nodes{i};
                end
            end
        end

        function statefulnames = getStatefulNodeNames(self)
            % STATEFULNAMES = GETSTATEFULNODENAMES()

            statefulnames = {};
            for i=1:self.getNumberOfNodes
                if self.nodes{i}.isStateful
                    statefulnames{end+1,1} = self.nodes{i}.name;
                end
            end
        end

        function M = getNumberOfNodes(self)
            % M = GETNUMBEROFNODES()

            M = length(self.nodes);
        end

        function S = getNumberOfStatefulNodes(self)
            % S = GETNUMBEROFSTATEFULNODES()

            S = sum(cellisa(self.nodes,'StatefulNode'));
        end

        function M = getNumberOfStations(self)
            % M = GETNUMBEROFSTATIONS()

            M = length(self.stations);
        end

        function R = getNumberOfClasses(self)
            % R = GETNUMBEROFCLASSES()

            R = length(self.classes);
        end

        function C = getNumberOfChains(self)
            % C = GETNUMBEROFCHAINS()

            sn = self.getStruct;
            C = sn.nchains;
        end

        function Dchain = getDemandsChain(self)
            % DCHAIN = GETDEMANDSCHAIN()
            sn_get_demands_chain(self.getStruct);
        end

        function self = setUsedLangFeature(self,className)
            % SELF = SETUSEDLANGFEATURE(SELF,CLASSNAME)

            self.usedFeatures.setTrue(className);
        end

        %% Add the components to the model
        addJobClass(self, customerClass);
        bool = addNode(self, node);
        fcr = addRegion(self, nodes);
        addLink(self, nodeA, nodeB);
        addLinks(self, nodeList);
        addItemSet(self, itemSet);

        node = getSource(self);
        node = getSink(self);

        function list = getStationIndexes(self)
            % LIST = GETSTATIONINDEXES()

            if self.hasStruct
                list = find(self.sn.isstation)';
            else
                % returns the ids of nodes that are stations
                list = find(cellisa(self.nodes, 'Station'))';
            end
        end

        function list = getIndexStatefulNodes(self)
            % LIST = GETINDEXSTATEFULNODES()

            % returns the ids of nodes that are stations
            list = find(cellisa(self.nodes, 'StatefulNode'))';
        end

        index = getIndexSourceStation(self);
        index = getIndexSourceNode(self);
        index = getIndexSinkNode(self);

        N = getNumberOfJobs(self);
        refstat = getReferenceStations(self);
        refclass = getReferenceClasses(self);
        sched = getStationScheduling(self);
        S = getStationServers(self);
        S = getStatefulServers(self);

        jsimwView(self)
        jsimgView(self)
        function view(self)
            % VIEW() Open the model in JSIMgraph
            jsimgView(self);
        end

        function modelView(self)
            % MODELVIEW() Open the model in ModelVisualizer
            jnetwork = JLINE.line_to_jline(self);
            jnetwork.plot();
        end

        function islld = isLimitedLoadDependent(self)
            islld = isempty(self.getStruct().lldscaling);
        end

        function [isvalid] = isStateValid(self)
            % [ISVALID] = ISSTATEVALID()

            isvalid = sn_is_state_valid(self.getStruct);
        end

        [state, priorStateSpace, stateSpace] = getState(self) % get initial state

        initFromAvgTableQLen(self, AvgTable)
        initFromAvgQLen(self, AvgQLen)
        initDefault(self, nodes)

        initFromMarginal(self, n, options) % n(i,r) : number of jobs of class r in node or station i (autodetected)
        initFromMarginalAndRunning(self, n, s, options) % n(i,r) : number of jobs of class r in node or station i (autodetected)
        initFromMarginalAndStarted(self, n, s, options) % n(i,r) : number of jobs of class r in node or station i (autodetected)

        [H,G] = getGraph(self)

        function mask = getClassSwitchingMask(self)
            % MASK = GETCLASSSWITCHINGMASK()

            mask = self.getStruct.csmask;
        end

        function printRoutingMatrix(self, onlyclass)
            % PRINTROUTINGMATRIX()
            if nargin==1
                sn_print_routing_matrix(self.getStruct);
            else
                sn_print_routing_matrix(self.getStruct, onlyclass);
            end
        end
    end

    % Private methods
    methods (Access = protected)

        function out = getModelNameExtension(self)
            % OUT = GETMODELNAMEEXTENSION()

            out = [getModelName(self), ['.', self.fileFormat]];
        end

        function self = initUsedFeatures(self)
            % SELF = INITUSEDFEATURES()

            % The list includes all classes but Model and Hidden or
            % Constant or Abstract or Solvers
            self.usedFeatures = SolverFeatureSet;
        end
    end

    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()

            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each handle
            for i=1:length(self.classes)
                clone.classes{i} = self.classes{i}.copy;
            end
            % Make a deep copy of each handle
            for i=1:length(self.nodes)
                clone.nodes{i} = self.nodes{i}.copy;
                if isa(clone.nodes{i},'Station')
                    clone.stations{i} = clone.nodes{i};
                end
                clone.connections = self.connections;
            end
        end
    end

    methods % wrappers
        function bool = hasFCFS(self)
            % BOOL = HASFCFS()

            bool = sn_has_fcfs(self.getStruct);
        end

        function bool = hasHomogeneousScheduling(self, strategy)
            % BOOL = HASHOMOGENEOUSSCHEDULING(STRATEGY)

            bool = sn_has_homogeneous_scheduling(self.getStruct, strategy);
        end

        function bool = hasDPS(self)
            % BOOL = HASDPS()

            bool = sn_has_dps(self.getStruct);
        end

        function bool = hasGPS(self)
            % BOOL = HASGPS()

            bool = sn_has_gps(self.getStruct);
        end

        function bool = hasINF(self)
            % BOOL = HASINF()

            bool = sn_has_inf(self.getStruct);
        end

        function bool = hasPS(self)
            % BOOL = HASPS()

            bool = sn_has_ps(self.getStruct);
        end

        function bool = hasSIRO(self)
            % BOOL = HASSIRO()

            bool = sn_has_siro(self.getStruct);
        end

        function bool = hasHOL(self)
            % BOOL = HASHOL()

            bool = sn_has_hol(self.getStruct);
        end

        function bool = hasLCFS(self)
            % BOOL = HASLCFS()

            bool = sn_has_lcfs(self.getStruct);
        end

        function bool = hasLCFSPR(self)
            % BOOL = HASLCFSPR()

            bool = sn_has_lcfs_pr(self.getStruct);
        end

        function bool = hasSEPT(self)
            % BOOL = HASSEPT()

            bool = sn_has_sept(self.getStruct);
        end

        function bool = hasLEPT(self)
            % BOOL = HASLEPT()

            bool = sn_has_lept(self.getStruct);
        end

        function bool = hasSJF(self)
            % BOOL = HASSJF()

            bool = sn_has_sjf(self.getStruct);
        end

        function bool = hasLJF(self)
            % BOOL = HASLJF()

            bool = sn_has_ljf(self.getStruct);
        end

        function bool = hasMultiClassFCFS(self)
            % BOOL = HASMULTICLASSFCFS()

            bool = sn_has_multi_class_fcfs(self.getStruct);
        end

        function bool = hasMultiClassHeterFCFS(self)
            % BOOL = HASMULTICLASSFCFS()

            bool = sn_has_multi_class_heter_fcfs(self.getStruct);
        end

        function bool = hasMultiServer(self)
            % BOOL = HASMULTISERVER()

            bool = sn_has_multi_server(self.getStruct);
        end

        function bool = hasSingleChain(self)
            % BOOL = HASSINGLECHAIN()

            bool = sn_has_single_chain(self.getStruct);
        end

        function bool = hasMultiChain(self)
            % BOOL = HASMULTICHAIN()

            bool = sn_has_multi_chain(self.getStruct);
        end

        function bool = hasSingleClass(self)
            % BOOL = HASSINGLECLASS()

            bool = sn_has_single_class(self.getStruct);
        end

        function bool = hasMultiClass(self)
            % BOOL = HASMULTICLASS()

            bool = sn_has_multi_class(self.getStruct);
        end
    end

    methods
        function bool = hasProductFormSolution(self)
            % BOOL = HASPRODUCTFORMSOLUTION()
            bool = sn_has_product_form(self.getStruct);
        end

        function plot(self)
            % PLOT()
            [~,H] = self.getGraph;
            H.Nodes.Name=strrep(H.Nodes.Name,'_','\_');
            h=plot(H,'EdgeLabel',H.Edges.Weight,'Layout','Layered');
            highlight(h,self.getNodeTypes==3,'NodeColor','r'); % class-switch nodes
        end

        function varargout = getMarkedCTMC( varargin )
            [varargout{1:nargout}] = getCTMC( varargin{:} );
        end

        function mctmc = getCTMC(self, par1, par2)
            if nargin<2
                options = SolverCTMC.defaultOptions;
            elseif ischar(par1)
                options = SolverCTMC.defaultOptions;
                options.(par1) = par2;
                options.cache = false;
            elseif isstruct(par1)
                options = par1;
            end
            solver = SolverCTMC(self,options);
            [infGen, eventFilt, ev] = solver.getInfGen();
            mctmc = MarkedMarkovProcess(infGen, eventFilt, ev, true);
            mctmc.setStateSpace(solver.getStateSpace);
        end
    end

    methods (Static)

        function model = tandemPs(lambda,D)
            % MODEL = TANDEMPS(LAMBDA,D)

            model = Network.tandemPsInf(lambda,D,[]);
        end

        function model = tandemPsInf(lambda,D,Z)
            % MODEL = TANDEMPSINF(LAMBDA,D,Z)

            if nargin<3%~exist('Z','var')
                Z = [];
            end
            M  = size(D,1);
            Mz = size(Z,1);
            strategy = {};
            for i=1:Mz
                strategy{i} = SchedStrategy.INF;
            end
            for i=1:M
                strategy{Mz+i} = SchedStrategy.PS;
            end
            model = Network.tandem(lambda,[Z;D],strategy);
        end

        function model = tandemFcfs(lambda,D)
            % MODEL = TANDEMFCFS(LAMBDA,D)

            model = Network.tandemFcfsInf(lambda,D,[]);
        end

        function model = tandemFcfsInf(lambda,D,Z)
            % MODEL = TANDEMFCFSINF(LAMBDA,D,Z)

            if nargin<3%~exist('Z','var')
                Z = [];
            end
            M  = size(D,1);
            Mz = size(Z,1);
            strategy = {};
            for i=1:Mz
                strategy{i} = SchedStrategy.INF;
            end
            for i=1:M
                strategy{Mz+i} = SchedStrategy.FCFS;
            end
            model = Network.tandem(lambda,[Z;D],strategy);
        end

        function model = tandem(lambda,D,strategy)
            % MODEL = TANDEM(LAMBDA,S,STRATEGY)

            % D(i,r) - mean service demand of class r at station i, equal
            % to mean service time since the topology is linear
            % lambda(r) - number of jobs of class r
            % station(i) - scheduling strategy at station i
            model = Network('Model');
            [M,R] = size(D);
            node{1} = Source(model, 'Source');
            for i=1:M
                switch SchedStrategy.toId(strategy{i})
                    case SchedStrategy.INF
                        node{end+1} = Delay(model, ['Station',num2str(i)]);
                    otherwise
                        node{end+1} = Queue(model, ['Station',num2str(i)], strategy{i});
                end
            end
            node{end+1} = Sink(model, 'Sink');
            P = cellzeros(R,R,M+2,M+2);
            for r=1:R
                jobclass{r} = OpenClass(model, ['Class',num2str(r)], 0);
                P{r,r} = circul(length(node)); P{r}(end,:) = 0;
            end
            for r=1:R
                node{1}.setArrival(jobclass{r}, Exp.fitMean(1/lambda(r)));
                for i=1:M
                    node{1+i}.setService(jobclass{r}, Exp.fitMean(D(i,r)));
                end
            end
            model.link(P);
        end

        function model = cyclicPs(N,D)
            % MODEL = CYCLICPS(N,D)
            M  = size(D,1);
            if nargin<4
                S = ones(M,1);
            end
            model = Network.cyclicPsInf(N,D,[],S);
        end

        function model = cyclicPsInf(N,D,Z,S)
            % MODEL = CYCLICPSINF(N,D,Z,S)
            if nargin<3
                Z = [];
            end
            M  = size(D,1);
            if nargin<4
                S = ones(M,1);
            end
            Mz = size(Z,1);
            strategy = cell(M+Mz,1);
            for i=1:Mz
                strategy{i} = SchedStrategy.INF;
            end
            for i=1:M
                strategy{Mz+i} = SchedStrategy.PS;
            end
            model = Network.cyclic(N,[Z;D],strategy,[Inf*ones(size(Z,1),1);S]);
        end

        function model = cyclicFcfs(N,D,S)
            % MODEL = CYCLICFCFS(N,D,S)
            M  = size(D,1);
            if nargin<4
                S = ones(M,1);
            end
            model = Network.cyclicFcfsInf(N,D,[],S);
        end

        function model = cyclicFcfsInf(N,D,Z,S)
            % MODEL = CYCLICFCFSINF(N,D,Z,S)

            if nargin<3%~exist('Z','var')
                Z = [];
            end
            M  = size(D,1);
            if nargin<4
                S = ones(M,1);
            end

            Mz = size(Z,1);
            strategy = {};
            for i=1:Mz
                strategy{i} = SchedStrategy.INF;
            end
            for i=1:M
                strategy{Mz+i} = SchedStrategy.FCFS;
            end
            model = Network.cyclic(N,[Z;D],strategy,[Inf*ones(size(Z,1),1);S]);
        end

        function model = cyclic(N,D,strategy,S)
            % MODEL = CYCLIC(N,D,STRATEGY,S)

            % L(i,r) - demand of class r at station i
            % N(r) - number of jobs of class r
            % strategy(i) - scheduling strategy at station i
            % S(i) - number of servers at station i
            model = Network('Model');
            [M,R] = size(D);
            node = cell(M,1);
            nQ = 0; nD = 0;
            for ist=1:M
                switch SchedStrategy.toId(strategy{ist})
                    case SchedStrategy.INF
                        nD = nD + 1;
                        node{ist} = Delay(model, ['Delay',num2str(nD)]);
                    otherwise
                        nQ = nQ + 1;
                        node{ist} = Queue(model, ['Queue',num2str(nQ)], strategy{ist});
                        node{ist}.setNumberOfServers(S(ist));
                end
            end
            P = cellzeros(R,M);
            jobclass = cell(R,1);
            for r=1:R
                jobclass{r} = ClosedClass(model, ['Class',num2str(r)], N(r), node{1}, 0);
                P{r,r} = circul(M);
            end
            for ist=1:M
                for r=1:R
                    node{ist}.setService(jobclass{r}, Exp.fitMean(D(ist,r)));
                end
            end
            model.link(P);
        end

        function P = serialRouting(varargin)
            % P = SERIALROUTING(VARARGIN)

            if length(varargin)==1
                varargin = varargin{1};
            end
            model = varargin{1}.model;
            P = zeros(model.getNumberOfNodes);
            for i=1:length(varargin)-1
                P(varargin{i},varargin{i+1})=1;
            end
            if ~isa(varargin{end},'Sink')
                P(varargin{end},varargin{1})=1;
            end
            P = P ./ repmat(sum(P,2),1,length(P));
            P(isnan(P)) = 0;
        end

        function printInfGen(Q,SS)
            % PRINTINFGEN(Q,SS)
            SolverCTMC.printInfGen(Q,SS);
        end

        function printEventFilt(sync,D,SS,myevents)
            % PRINTEVENTFILT(SYNC,D,SS,MYEVENTS)
            SolverCTMC.printEventFilt(sync,D,SS,myevents);
        end

    end
end
