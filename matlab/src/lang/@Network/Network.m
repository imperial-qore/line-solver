classdef Network < Model
    % An extended queueing network model.
    %
    % Copyright (c) 2012-2024, Imperial College London
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
    end

    methods % link, reset, refresh methods
        self = link(self, P)
        self = linkFromNodeRoutingMatrix(self, Pnodes);
        [loggerBefore,loggerAfter] = linkAndLog(self, nodes, classes, P, wantLogger, logPath)
        [loggerBefore,loggerAfter] = linkNetworkAndLog(self, nodes, classes, P, wantLogger, logPath)% obsolete - old name
        sanitize(self);
        self = removeClass(self, jobclass);

        reset(self, resetState)
        resetHandles(self)
        resetModel(self, resetState)
        nodes = resetNetwork(self, deleteCSnodes)
        resetStruct(self)

        refreshStruct(self, hard);
        [rates, scv, hasRateChanged, hasSCVChanged] = refreshRates(self, statSet, classSet);
        [ph, mu, phi, phases] = refreshServicePhases(self, statSet, classSet);
        proctypes = refreshServiceTypes(self);
        [rt, rtfun, rtnodes] = refreshRoutingMatrix(self, rates);
        [lt] = refreshLST(self, statSet, classSet);
        sync = refreshSync(self);
        classprio = refreshPriorities(self);
        [sched, schedid, schedparam] = refreshScheduling(self);
        [rates, mu, phi, phases] = refreshArrival(self);
        [rates, scv, mu, phi, phases] = refreshService(self, statSet, classSet);
        [chains, visits, rt] = refreshChains(self, propagate)
        [cap, classcap] = refreshCapacity(self);
        nvars = refreshLocalVars(self);
        [nonfjmodel, fjclassmap, fjforkmap, fanout] = approxForkJoins(self, forkLambda);
    end

    % PUBLIC METHODS
    methods (Access=public)

        %Constructor
        function self = Network(modelName)
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

        function index = getIndexOpenClasses(self)
            % INDEX = GETINDEXOPENCLASSES()

            index = find(isinf(getNumberOfJobs(self)))';
        end

        function index = getIndexClosedClasses(self)
            % INDEX = GETINDEXCLOSEDCLASSES()

            index = find(isfinite(getNumberOfJobs(self)))';
        end

        chain = getClassChain(self, className)
        c = getClassChainIndex(self, className)
        classNames = getClassNames(self)

        nodeNames = getNodeNames(self)
        nodeTypes = getNodeTypes(self)

        P = initRoutingMatrix(self)
        P = allCyclicRoutingMatrix(self)

        rtTypes = getRoutingStrategies(self)
        ind = getNodeIndex(self, name)
        lldScaling = getLimitedLoadDependence(self)
        lcdScaling = getLimitedClassDependence(self)
        simodel = getStateIndependentModel(self);

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
                for i=self.getIndexStations
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

        function P=summaryChains(self)
            if ~self.hasStruct()
                self.getStruct;
            end
            G = digraph(self.sn.rtnodes);
            for c=1:self.sn.nchains
                inchain = find(self.sn.chains(c,:));
                if self.sn.refclass(c)>0
                    root = (self.sn.stationToNode(self.sn.refstat(inchain(1)))-1)+self.sn.refclass(c);
                    [~,E]=G.dfsearch(root,'edgetonew');
                    [~,E2]=G.dfsearch(root,'edgetodiscovered');
                    [~,E3]=G.dfsearch(root,'edgetofinished');
                    E(end+1:end+size(E2,1),:)=E2;
                    E(end+1:end+size(E3,1),:)=E3;
                    for j=1:size(E,1)                        
                        rfrom=mod(G.Edges.EndNodes(E(j),1),self.sn.nclasses);
                        if rfrom==0
                            rfrom = self.sn.nclasses;
                        end
                        ifrom=(G.Edges.EndNodes(E(j),1)-mod(G.Edges.EndNodes(E(j),1),self.sn.nclasses))/self.sn.nclasses + 1;
                        rto=mod(G.Edges.EndNodes(E(j),2),self.sn.nclasses);
                        if rto==0
                            rto = self.sn.nclasses;
                        end
                        ito=(G.Edges.EndNodes(E(j),2)-mod(G.Edges.EndNodes(E(j),2),self.sn.nclasses))/self.sn.nclasses + 1;
                        P{c}{j,1} = self.sn.classnames{rfrom};
                        P{c}{j,2} = self.sn.nodenames{ifrom};
                        P{c}{j,3} = self.sn.classnames{rto};
                        P{c}{j,4} = self.sn.nodenames{ito};
                    end
                    P{c} = sortrows(P{c},'ascend');
                else
                    line_warning(mfilename,'Chain %d has no reference class set, skipping.\n',c);
                end
            end
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
        end

        function [D,Z] = getDemands(self)
            % [D,Z]= GETDEMANDS()
            %
            % Outputs:
            % D: service demands
            % Z: think times

            [~,D,~,Z,~,~] = snGetProductFormParams(self.getStruct);
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

            [lambda,D,N,Z,mu,S] = snGetProductFormParams(self.getStruct);
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

            [lambda,D,N,Z,mu,S]=  snGetProductFormChainParams(self.getStruct);
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
            snGetDemandsChain(self.getStruct);
        end

        % setUsedFeatures : records that a certain language feature has been used
        function self = setUsedFeatures(self,className)
            % SELF = SETUSEDFEATURES(SELF,CLASSNAME)

            self.usedFeatures.setTrue(className);
        end

        %% Add the components to the model
        addJobClass(self, customerClass);
        bool = addNode(self, node);
        fcr = addRegion(self, nodes);
        addLink(self, nodeA, nodeB);
        addLinks(self, nodeList);
        addItemSet(self, itemSet);
        self = disableMetric(self, Y);
        self = enableMetric(self, Y);

        node = getSource(self);
        node = getSink(self);

        function list = getDummys(self)
            % LIST = GETDUMMYS()

            list = find(cellisa(self.nodes, 'Passage'))';
        end

        function list = getIndexStations(self)
            % LIST = GETINDEXSTATIONS()

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

        jsimwView(self)
        jsimgView(self)

        function [ni, nir, sir, kir] = initToMarginal(self)
            % [NI, NIR, SIR, KIR] = INITTOMARGINAL()

            [ni, nir, sir, kir] = snInitToMarginal(self.getStruct);
        end

        function islld = isLimitedLoadDependent(self)
            sn = self.getStruct;
            if isempty(sn.lldscaling)
                islld = false;
            else
                islld = true;
            end
        end

        function [isvalid] = isStateValid(self)
            % [ISVALID] = ISSTATEVALID()

            isvalid = snIsStateValid(self.getStruct);
        end

        [initialStateAggr] = getStateAggr(self) % get initial state
        [initialState, priorInitialState] = getState(self) % get initial state

        initFromAvgTableQLen(self, AvgTable)
        initFromAvgQLen(self, AvgQLen)
        initDefault(self, nodes)
        initFromMarginal(self, n, options) % n(i,r) : number of jobs of class r in node i
        initFromMarginalAndRunning(self, n, s, options) % n(i,r) : number of jobs of class r in node i
        initFromMarginalAndStarted(self, n, s, options) % n(i,r) : number of jobs of class r in node i

        [H,G] = getGraph(self)

        function mask = getClassSwitchingMask(self)
            % MASK = GETCLASSSWITCHINGMASK()

            mask = self.getStruct.csmask;
        end

        function printRoutingMatrix(self, onlyclass)
            % PRINTROUTINGMATRIX()
            if nargin==1
                snPrintRoutingMatrix(self.getStruct);
            else
                snPrintRoutingMatrix(self.getStruct, onlyclass);
            end
        end

        function [taggedModel, taggedJob] = tagChain(self, chain, jobclass, suffix)
            % the tagged job will be removed from the initial
            % population of JOBCLASS
            if nargin<4 || isempty(suffix)
                suffix = '.tagged';
            end
            if nargin<3
                jobclass = chain.classes{1};
            end
            I = self.getNumberOfNodes;
            R = self.getNumberOfClasses;
            taggedModel = self.copy;

            % we don't use rtNodesByClass because it contains the
            % fictitious class switching nodes
            Plinked = taggedModel.getLinkedRoutingMatrix;
            if ~iscell(Plinked) || isempty(Plinked)
                line_error(mfilename, 'getCdfRespT requires the original model to be linked with a routing matrix defined as a cell array P{r,s} for every class pair (r,s).');
            end

            taggedModel.resetNetwork; % resets cs Nodes as well
            taggedModel.reset(true);

            chainIndexes = cell2mat(chain.index);
            for r=chainIndexes
                % create a tagged class
                taggedModel.classes{end+1,1} = taggedModel.classes{r}.copy;
                taggedModel.classes{end,1}.index=length(taggedModel.classes);
                taggedModel.classes{end,1}.name=[taggedModel.classes{r,1}.name,suffix];
                if r==jobclass.index
                    taggedModel.classes{end}.population = 1;
                else
                    taggedModel.classes{end}.population = 0;
                end

                % clone station sections for tagged class
                for m=1:length(taggedModel.nodes)
                    taggedModel.stations{m}.output.outputStrategy{end+1} = taggedModel.stations{m}.output.outputStrategy{r};
                end

                for m=1:length(taggedModel.stations)
                    if self.stations{m}.server.serviceProcess{r}{end}.isDisabled
                        taggedModel.stations{m}.input.inputJobClasses(1,end+1) = {[]};
                        taggedModel.stations{m}.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{end}{end}.copy;
                        taggedModel.stations{m}.server.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r};
                        taggedModel.stations{m}.server.serviceProcess{end}{end}=taggedModel.stations{m}.server.serviceProcess{r}{end}.copy;
                        taggedModel.stations{m}.schedStrategyPar(end+1) = 0;
                        taggedModel.stations{m}.dropRule(1,end+1) = -1;
                        taggedModel.stations{m}.classCap(1,r) = 0;
                        taggedModel.stations{m}.classCap(1,end+1) = 0;
                    else
                        taggedModel.stations{m}.input.inputJobClasses(1,end+1) = {taggedModel.stations{m}.input.inputJobClasses{1,r}};
                        taggedModel.stations{m}.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r}{end}.copy;
                        taggedModel.stations{m}.server.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r};
                        taggedModel.stations{m}.server.serviceProcess{end}{end}=taggedModel.stations{m}.server.serviceProcess{r}{end}.copy;
                        taggedModel.stations{m}.schedStrategyPar(end+1) = taggedModel.stations{m}.schedStrategyPar(r);
                        taggedModel.stations{m}.classCap(1,end+1) = 1;
                        taggedModel.stations{m}.dropRule(1,end+1) = -1;
                        taggedModel.stations{m}.classCap(1,r) = taggedModel.stations{m}.classCap(r) - 1;
                    end
                end
            end

            taggedModel.classes{jobclass.index,1}.population = taggedModel.classes{jobclass.index}.population - 1;

            for ir=1:length(chainIndexes)
                r = chainIndexes(ir);
                for is=1:length(chainIndexes)
                    s = chainIndexes(is);
                    Plinked{R+ir,R+is} = Plinked{r,s};
                end
            end
            Rp = taggedModel.getNumberOfClasses;
            for r=1:Rp
                for s=1:Rp
                    if isempty(Plinked{r,s})
                        Plinked{r,s} = zeros(I);
                    end
                end
            end
            taggedModel.sn = [];
            taggedModel.link(Plinked);
            taggedModel.reset(true);
            taggedModel.refreshStruct(true);
            taggedModel.initDefault;
            tchains = taggedModel.getChains;
            taggedJob = tchains{end};
        end

    end

    % Private methods
    methods (Access = 'private')

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

            bool = snHasFCFS(self.getStruct);
        end

        function bool = hasHomogeneousScheduling(self, strategy)
            % BOOL = HASHOMOGENEOUSSCHEDULING(STRATEGY)

            bool = snHasHomogeneousScheduling(self.getStruct, strategy);
        end

        function bool = hasDPS(self)
            % BOOL = HASDPS()

            bool = snHasDPS(self.getStruct);
        end

        function bool = hasGPS(self)
            % BOOL = HASGPS()

            bool = snHasGPS(self.getStruct);
        end

        function bool = hasINF(self)
            % BOOL = HASINF()

            bool = snHasINF(self.getStruct);
        end

        function bool = hasPS(self)
            % BOOL = HASPS()

            bool = snHasPS(self.getStruct);
        end

        function bool = hasSIRO(self)
            % BOOL = HASSIRO()

            bool = snHasSIRO(self.getStruct);
        end

        function bool = hasHOL(self)
            % BOOL = HASHOL()

            bool = snHasHOL(self.getStruct);
        end

        function bool = hasLCFS(self)
            % BOOL = HASLCFS()

            bool = snHasLCFS(self.getStruct);
        end

        function bool = hasLCFSPR(self)
            % BOOL = HASLCFSPR()

            bool = snHasLCFSPR(self.getStruct);
        end
        
        function bool = hasSEPT(self)
            % BOOL = HASSEPT()

            bool = snHasSEPT(self.getStruct);
        end

        function bool = hasLEPT(self)
            % BOOL = HASLEPT()

            bool = snHasLEPT(self.getStruct);
        end

        function bool = hasSJF(self)
            % BOOL = HASSJF()

            bool = snHasSJF(self.getStruct);
        end

        function bool = hasLJF(self)
            % BOOL = HASLJF()

            bool = snHasLJF(self.getStruct);
        end

        function bool = hasMultiClassFCFS(self)
            % BOOL = HASMULTICLASSFCFS()

            bool = snHasMultiClassFCFS(self.getStruct);
        end

        function bool = hasMultiClassHeterFCFS(self)
            % BOOL = HASMULTICLASSFCFS()

            bool = snHasMultiClassHeterFCFS(self.getStruct);
        end

        function bool = hasMultiServer(self)
            % BOOL = HASMULTISERVER()

            bool = snHasMultiServer(self.getStruct);
        end

        function bool = hasSingleChain(self)
            % BOOL = HASSINGLECHAIN()

            bool = snHasSingleChain(self.getStruct);
        end

        function bool = hasMultiChain(self)
            % BOOL = HASMULTICHAIN()

            bool = snHasMultiChain(self.getStruct);
        end

        function bool = hasSingleClass(self)
            % BOOL = HASSINGLECLASS()

            bool = snHasSingleClass(self.getStruct);
        end

        function bool = hasMultiClass(self)
            % BOOL = HASMULTICLASS()

            bool = snHasMultiClass(self.getStruct);
        end

        function print(self)
            LINE2SCRIPT(self)
        end
    end

    methods
        function bool = hasProductFormSolution(self)
            % BOOL = HASPRODUCTFORMSOLUTION()
            bool = snHasProductForm(self.getStruct);
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
            mctmc = MarkedCTMC(infGen, eventFilt, ev, true);
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
                    case SchedStrategy.ID_INF
                        node{end+1} = DelayStation(model, ['Station',num2str(i)]);
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

        function model = productForm(N,D,Z)
            if nargin<3
                Z = [];
            end
            model = Network.cyclicPsInf(N,D,Z);
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
            for i=1:M
                switch SchedStrategy.toId(strategy{i})
                    case SchedStrategy.ID_INF
                        nD = nD + 1;
                        node{i} = DelayStation(model, ['Delay',num2str(nD)]);
                    otherwise
                        nQ = nQ + 1;
                        node{i} = Queue(model, ['Queue',num2str(nQ)], strategy{i});
                        node{i}.setNumberOfServers(S(i));
                end
            end
            P = cellzeros(R,M);
            jobclass = cell(R,1);
            for r=1:R
                jobclass{r} = ClosedClass(model, ['Class',num2str(r)], N(r), node{1}, 0);
                P{r,r} = circul(M);
            end
            for i=1:M
                for r=1:R
                    node{i}.setService(jobclass{r}, Exp.fitMean(D(i,r)));
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
