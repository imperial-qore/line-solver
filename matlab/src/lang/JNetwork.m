classdef JNetwork < Model
    % JLINE extended queueing network model.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Access=public)
        obj; % java object
    end

    % PUBLIC METHODS
    methods (Access=public)

        %Constructor
        function self = JNetwork(name)
            % SELF = NETWORK(MODELNAME)
            self@Model(name); % model is the model's name
            self.obj = jline.lang.Network(name);
        end

        function sn = getStruct(self, wantInitialState) % get abritrary representation
            if nargin<2
                wantInitialState = false;
            end
            sn = JLINE.from_jline_struct(self.obj, self.obj.getStruct(wantInitialState));
        end

        function R = getNumberOfClasses(self)
            % R = GETNUMBEROFCLASSES()

            R = self.obj.getNumberOfClasses();
        end

        function M = getNumberOfStations(self)
            % M = GETNUMBEROFSTATIONS()

            M = self.obj.getNumberOfStations();
        end


        function I = getNumberOfNodes(self)
            % I = GETNUMBEROFNODES()

            I = self.obj.getNumberOfNodes();
        end

        function bool = hasFork(self)
            % to be changed
            bool = false;
        end

        function ind = getNodeIndex(self, name)
            if ischar(name)
                ind = self.obj.getNodeByName(name);
            else
                ind = self.obj.getNodeIndex(name.obj);
            end
        end

        function node = getNodeByName(self, name)
            % NODE = GETNODEBYNAME(NAME)
            
            javaNode = self.obj.getNodeByName(name);
            if isempty(javaNode)
                node = NaN;
            else
                % Create a wrapper MATLAB object for the Java node
                % We need to determine the node type and create appropriate wrapper
                nodeClassName = char(javaNode.getClass().getSimpleName());
                
                % Create the MATLAB wrapper object, then override the Java object
                switch nodeClassName
                    case 'Queue'
                        node = Queue(self, char(javaNode.getName()), SchedStrategy.FCFS);  % Default schedule strategy
                    case 'Source'
                        node = Source(self, char(javaNode.getName()));
                    case 'Sink'
                        node = Sink(self, char(javaNode.getName()));
                    case 'Router'
                        node = Router(self, char(javaNode.getName()));
                    case 'Fork'
                        node = Fork(self, char(javaNode.getName()));
                    case 'Join'
                        node = Join(self, char(javaNode.getName()));
                    case 'Cache'
                        % Cache requires additional parameters, using defaults
                        node = Cache(self, char(javaNode.getName()), 1, [1], ReplacementStrategy.LRU);
                    case 'Logger'
                        node = Logger(self, char(javaNode.getName()), 'default.log');
                    case 'ClassSwitch'
                        node = ClassSwitch(self, char(javaNode.getName()));
                    case 'Place'
                        node = Place(self, char(javaNode.getName()));
                    case 'Transition'
                        node = Transition(self, char(javaNode.getName()));
                    otherwise
                        % Fallback - create generic Node wrapper
                        node = Node(char(javaNode.getName()));
                        node.setModel(self);
                end
                
                % Override the Java object with the existing one from the network
                node.obj = javaNode;
                node.index = self.obj.getNodeIndex(javaNode);
            end
        end

        function self = link(self, P)
            if isa(P,'jline.lang.RoutingMatrix')
                self.obj.link(P);
            else
                nodes = self.obj.getNodes();
                classes = self.obj.getClasses();
                routing_matrix = jline.lang.RoutingMatrix(self.obj, classes, nodes);
                I = self.obj.getNumberOfNodes;
                R = self.obj.getNumberOfClasses;

                if size(P,1) == 1
                    for i = 1:I
                        for j = 1:I
                            for r=1:R
                                for s=1:R
                                    if P{r,s}(i,j) > 0
                                        routing_matrix.addConnection(nodes.get(i-1), nodes.get(j-1), classes.get(r-1), classes.get(s-1), P{r,s}(i,j));
                                    end
                                end
                            end
                        end
                    end
                else % if a single matrix
                    for i = 1:I
                        for j = 1:I
                            if P(i,j) > 0
                                routing_matrix.addConnection(nodes.get(i-1), nodes.get(j-1), classes.get(0), P(i,j));
                            end
                        end
                    end
                end
                self.obj.link(routing_matrix);
            end
        end

        function P = initRoutingMatrix(self)
            M = self.getNumberOfNodes();
            R = self.getNumberOfClasses();
            P = cellzeros(R,R,M,M);
            P = RoutingMatrix(P);
        end

        function self = setChecks(self, bool)
            self.obj.setChecks(bool);
        end

        function addLinks(self, nodesList)
            % ADDLINKS(NODESLIST)

            % Copyright (c) 2012-2026, Imperial College London
            % All rights reserved.
            for i=1:size(nodesList,1)
                self.obj.addLink(self.obj.getNodeByIndex(nodesList(i,1)-1), self.obj.getNodeByIndex(nodesList(i,2)-1));
            end
        end

        function addLink(self, node1, node2)
            % ADDLINK(NODESLIST)

            % Copyright (c) 2012-2026, Imperial College London
            % All rights reserved.
            self.obj.addLink(node1.obj, node2.obj);
        end

        function reset(self)
            self.obj.reset();
        end

        function jsimgView(self)
            self.obj.jsimgView();
        end

        % Additional methods for coherence with MNetwork
        
        function nodes = getNodes(self)
            % NODES = GETNODES()
            javaNodes = self.obj.getNodes();
            nodes = {};
            for i = 0:javaNodes.size()-1
                jNode = javaNodes.get(i);
                nodes{end+1} = self.getNodeByName(char(jNode.getName()));
            end
        end
        
        function classes = getClasses(self)
            % CLASSES = GETCLASSES()
            javaClasses = self.obj.getClasses();
            classes = {};
            for i = 0:javaClasses.size()-1
                jClass = javaClasses.get(i);
                % Note: Java object wrapping for classes not implemented
                classes{end+1} = jClass;
            end
        end
        
        function classNames = getClassNames(self)
            % CLASSNAMES = GETCLASSNAMES()
            javaClasses = self.obj.getClasses();
            classNames = {};
            for i = 0:javaClasses.size()-1
                classNames{end+1} = char(javaClasses.get(i).getName());
            end
        end
        
        function nodeNames = getNodeNames(self)
            % NODENAMES = GETNODENAMES()
            javaNodes = self.obj.getNodes();
            nodeNames = {};
            for i = 0:javaNodes.size()-1
                nodeNames{end+1} = char(javaNodes.get(i).getName());
            end
        end
        
        function nodeTypes = getNodeTypes(self)
            % NODETYPES = GETNODETYPES()
            javaNodeTypes = self.obj.getNodeTypes();
            nodeTypes = {};
            for i = 0:javaNodeTypes.size()-1
                nodeTypes{end+1} = char(javaNodeTypes.get(i).toString());
            end
        end
        
        function stationNames = getStationNames(self)
            % STATIONNAMES = GETSTATIONNAMES()
            % Note: Java object does not support getStationNames method directly
            nodes = self.getNodes();
            stationNames = {};
            for i = 1:length(nodes)
                if isa(nodes{i}, 'Station')
                    stationNames{end+1} = nodes{i}.getName();
                end
            end
        end
        
        function stationIndex = getStationIndex(self, name)
            % STATIONINDEX = GETSTATIONINDEX(NAME)
            if isa(name,'Node')
                node = name;
                name = node.getName();
            end
            stationIndex = find(cellfun(@(c) strcmp(c,name),self.getStationNames));
        end
        
        function classIndex = getClassIndex(self, name)
            % CLASSINDEX = GETCLASSINDEX(NAME)
            if isa(name,'JobClass')
                jobclass = name;
                name = jobclass.getName();
            end
            classIndex = find(cellfun(@(c) strcmp(c,name),self.getClassNames));
        end
        
        function station = getStationByName(self, name)
            % STATION = GETSTATIONBYNAME(NAME)
            node = self.getNodeByName(name);
            if isa(node, 'Station')
                station = node;
            else
                station = NaN;
            end
        end
        
        function class = getClassByName(self, name)
            % CLASS = GETCLASSBYNAME(NAME)
            javaClass = self.obj.getClassByName(name);
            if isempty(javaClass)
                class = NaN;
            else
                % Note: Java object wrapping for classes not implemented
                class = javaClass;
            end
        end
        
        function node = getNodeByIndex(self, idx)
            % NODE = GETNODEBYINDEX(IDX)
            javaNode = self.obj.getNodeByIndex(idx-1); % Java uses 0-based indexing
            if isempty(javaNode)
                node = NaN;
            else
                node = self.getNodeByName(char(javaNode.getName()));
            end
        end
        
        function station = getStationByIndex(self, idx)
            % STATION = GETSTATIONBYINDEX(IDX)
            stations = self.getStationNames();
            if idx > 0 && idx <= length(stations)
                station = self.getStationByName(stations{idx});
            else
                station = NaN;
            end
        end
        
        function class = getClassByIndex(self, idx)
            % CLASS = GETCLASSBYINDEX(IDX)
            javaClass = self.obj.getClassByIndex(idx-1); % Java uses 0-based indexing
            if isempty(javaClass)
                class = NaN;
            else
                % Note: Java object wrapping for classes not implemented
                class = javaClass;
            end
        end
        
        function C = getNumberOfChains(self)
            % C = GETNUMBEROFCHAINS()
            sn = self.getStruct();
            C = sn.nchains;
        end
        
        function S = getNumberOfStatefulNodes(self)
            % S = GETNUMBEROFSTATEFULNODES()
            % Note: Java object does not support this method directly
            nodes = self.getNodes();
            S = sum(cellfun(@(n) isa(n, 'StatefulNode'), nodes));
        end
        
        function list = getStationIndexes(self)
            % LIST = GETSTATIONINDEXES()
            nodes = self.getNodes();
            list = find(cellfun(@(n) isa(n, 'Station'), nodes))';
        end
        
        function list = getIndexStatefulNodes(self)
            % LIST = GETINDEXSTATEFULNODES()
            nodes = self.getNodes();
            list = find(cellfun(@(n) isa(n, 'StatefulNode'), nodes))';
        end
        
        function index = getIndexSourceStation(self)
            % INDEX = GETINDEXSOURCESTATION()
            % Note: Java object does not support this method directly
            nodes = self.getNodes();
            for i = 1:length(nodes)
                if isa(nodes{i}, 'Source')
                    index = i;
                    return;
                end
            end
            index = [];
        end
        
        function index = getIndexSourceNode(self)
            % INDEX = GETINDEXSOURCENODE()
            % Note: Java object does not support this method directly
            nodes = self.getNodes();
            for i = 1:length(nodes)
                if isa(nodes{i}, 'Source')
                    index = i;
                    return;
                end
            end
            index = [];
        end
        
        function index = getIndexSinkNode(self)
            % INDEX = GETINDEXSINKNODE()
            % Note: Java object does not support this method directly
            nodes = self.getNodes();
            for i = 1:length(nodes)
                if isa(nodes{i}, 'Sink')
                    index = i;
                    return;
                end
            end
            index = [];
        end
        
        function node = getSource(self)
            % NODE = GETSOURCE()
            idx = self.getIndexSourceNode();
            if ~isempty(idx)
                node = self.getNodeByIndex(idx);
            else
                node = NaN;
            end
        end
        
        function node = getSink(self)
            % NODE = GETSINK()
            idx = self.getIndexSinkNode();
            if ~isempty(idx)
                node = self.getNodeByIndex(idx);
            else
                node = NaN;
            end
        end
        
        function N = getNumberOfJobs(self)
            % N = GETNUMBEROFJOBS()
            javaClasses = self.obj.getClasses();
            N = zeros(1, javaClasses.size());
            for i = 0:javaClasses.size()-1
                jClass = javaClasses.get(i);
                if strcmp(char(jClass.getClass().getSimpleName()), 'ClosedClass')
                    N(i+1) = jClass.getPopulation();
                else
                    N(i+1) = Inf;
                end
            end
        end
        
        function bool = hasOpenClasses(self)
            % BOOL = HASOPENCLASSES()
            N = self.getNumberOfJobs();
            bool = any(isinf(N));
        end
        
        function bool = hasClosedClasses(self)
            % BOOL = HASCLOSEDCLASSES()
            N = self.getNumberOfJobs();
            bool = any(isfinite(N));
        end

        function bool = isMatlabNative(self)
            % BOOL = ISMATLABNATIVE()
            %
            % Returns false for Java (JNetwork) implementation
            
            bool = false;
        end

        function bool = isJavaNative(self)
            % BOOL = ISJAVANATIVE()
            %
            % Returns true for Java (JNetwork) implementation
            
            bool = true;
        end
        
        function index = getIndexOpenClasses(self)
            % INDEX = GETINDEXOPENCLASSES()
            N = self.getNumberOfJobs();
            index = find(isinf(N))';
        end
        
        function index = getIndexClosedClasses(self)
            % INDEX = GETINDEXCLOSEDCLASSES()
            N = self.getNumberOfJobs();
            index = find(isfinite(N))';
        end
        
        function bool = hasClassSwitching(self)
            % BOOL = HASCLASSSWITCHING()
            nodes = self.getNodes();
            bool = any(cellfun(@(n) isa(n,'ClassSwitch'), nodes));
        end
        
        function bool = hasJoin(self)
            % BOOL = HASJOIN()
            nodes = self.getNodes();
            bool = any(cellfun(@(n) isa(n,'Join'), nodes));
        end
        
        function [M,R] = getSize(self)
            % [M,R] = GETSIZE()
            M = self.getNumberOfNodes();
            R = self.getNumberOfClasses();
        end
        
        function used = getUsedLangFeatures(self)
            % USED = GETUSEDLANGFEATURES()
            javaFeatureSet = self.obj.getUsedLangFeatures();
            
            % Convert Java FeatureSet to MATLAB struct
            % Note: This is a simplified conversion - may need refinement
            % based on the actual FeatureSet structure
            used = struct();
            
            % If the Java FeatureSet has accessible fields/methods,
            % we would convert them here. For now, return the Java object
            % This may need custom conversion logic depending on FeatureSet implementation
            try
                % Attempt to get feature names if FeatureSet has such methods
                % This is a placeholder - actual implementation depends on FeatureSet API
                used = javaFeatureSet;
            catch
                % Fallback - create empty struct
                used = struct();
            end
        end
        
        function summary(self)
            % SUMMARY()
            self.obj.summary();
        end
        
        function [D,Z] = getDemands(self)
            % [D,Z] = GETDEMANDS()
            ret = self.obj.getDemands();
            D = JLINE.from_jline_matrix(ret.D);
            Z = JLINE.from_jline_matrix(ret.Z);
        end
        
        function [lambda,D,N,Z,mu,S] = getProductFormParameters(self)
            % [LAMBDA,D,N,Z,MU,S] = GETPRODUCTFORMPARAMETERS()
            ret = self.obj.getProductFormParameters();
            lambda = JLINE.from_jline_matrix(ret.lambda);
            D = JLINE.from_jline_matrix(ret.D);
            N = JLINE.from_jline_matrix(ret.N);
            Z = JLINE.from_jline_matrix(ret.Z);
            mu = JLINE.from_jline_matrix(ret.mu);
            S = JLINE.from_jline_matrix(ret.S);
        end
        
        function printRoutingMatrix(self, onlyclass)
            % PRINTROUTINGMATRIX(ONLYCLASS)
            self.obj.printRoutingMatrix();
        end
        
        function [rt,rtNodes,connections,chains,rtNodesByClass,rtNodesByStation] = getRoutingMatrix(self, arvRates)
            % [RT,RTNODES,CONNECTIONS,CHAINS,RTNODEBYCLASS,RTNODEBYSTATION] = GETROUTINGMATRIX(ARVRATES)
            % Note: Java object does not support this method directly
            error('getRoutingMatrix method not supported by Java object');
        end
        
        function P = getLinkedRoutingMatrix(self)
            % P = GETLINKEDROUTINGMATRIX()
            javaP = self.obj.getLinkedRoutingMatrix();
            % Convert Java Map<JobClass, Map<JobClass, Matrix>> to MATLAB cell array
            P = {};
            % Note: This is a complex conversion - simplified implementation
            % May need custom conversion logic based on how MATLAB expects this data
            P = javaP; % Direct return for now - may need refinement
        end
        
        function connections = getConnectionMatrix(self)
            % CONNECTIONS = GETCONNECTIONMATRIX()
            javaConnections = self.obj.getConnectionMatrix();
            connections = JLINE.from_jline_matrix(javaConnections);
        end
        
        function mask = getClassSwitchingMask(self)
            % MASK = GETCLASSSWITCHINGMASK()
            javaMask = self.obj.getClassSwitchingMask();
            mask = JLINE.from_jline_matrix(javaMask);
        end
        
        function bool = hasProductFormSolution(self)
            % BOOL = HASPRODUCTFORMSOLUTION()
            bool = self.obj.hasProductFormSolution();
        end
        
        function plot(self)
            % PLOT() - Display network as TikZ diagram in PDF viewer
            % Requires pdflatex to be installed on the system
            self.obj.plot();
        end
        
        function jsimwView(self)
            % JSIMWVIEW()
            self.obj.jsimwView();
        end
        
        function view(self)
            % VIEW() Open the model in JSIMgraph
            self.jsimgView();
        end

        function modelView(self)
            % MODELVIEW() Open the model in ModelVisualizer
            self.obj.plot();
        end
        
        % Save methods (for consistency with MNetwork/SolverJMT)
        function saveAsJSIM(self, filename)
            % SAVEASJSIM(FILENAME)
            % Note: Java object does not support direct JSIM export
            error('saveAsJSIM method not supported by Java object. Use SolverJMT for JSIM export.');
        end
        
        function saveAsJMVA(self, filename)
            % SAVEASJMVA(FILENAME)
            % Note: Java object does not support direct JMVA export
            error('saveAsJMVA method not supported by Java object. Use SolverJMT for JMVA export.');
        end

    end

    methods (Static)


        function model = tandemPs(lambda, D, S)
            lambda = JLINE.from_line_matrix(lambda);
            D = JLINE.from_line_matrix(D);
            S = JLINE.from_line_matrix(S);
            model = jline.lang.Network.tandemPs(lambda, D, S);
            model = JNetwork(model);
        end

        function model = tandemPsInf(varargin)
            lambda = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            if nargin == 2
                model = jline.lang.Network.tandemPsInf(lambda, D);
            elseif nargin == 3
                Z = JLINE.from_line_matrix(varargin{3});
                model = jline.lang.Network.tandemPsInf(lambda, D, Z);
            elseif nargin == 4
                Z = JLINE.from_line_matrix(varargin{3});
                S = JLINE.from_line_matrix(varargin{4});
                model = jline.lang.Network.tandemPsInf(lambda, D, Z, S);
            else
                error('Invalid number of arguments for tandemPsInf.');
            end
            model = JNetwork(model);
        end

        function model = tandemFcfs(varargin)
            lambda = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            S = JLINE.from_line_matrix(varargin{3});
            model = jline.lang.Network.tandemFcfs(lambda, D, S);
            model = JNetwork(model);
        end

        function model = tandemFcfsInf(varargin)
            lambda = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            if nargin == 2
                model = jline.lang.Network.tandemFcfsInf(lambda, D);
            elseif nargin == 3
                Z = JLINE.from_line_matrix(varargin{3});
                model = jline.lang.Network.tandemFcfsInf(lambda, D, Z);
            elseif nargin == 4
                Z = JLINE.from_line_matrix(varargin{3});
                S = JLINE.from_line_matrix(varargin{4});
                model = jline.lang.Network.tandemFcfsInf(lambda, D, Z, S);
            else
                error('Invalid number of arguments for tandemFcfsInf.');
            end
            model = JNetwork(model);
        end

        function model = tandem(lambda, D, strategy, S)
            lambda = JLINE.from_line_matrix(lambda);
            D = JLINE.from_line_matrix(D);
            S = JLINE.from_line_matrix(S);
            model = jline.lang.Network.tandem(lambda, D, strategy, S);
            model = JNetwork(model);
        end

        function model = cyclicPs(varargin)
            N = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            if nargin == 2
                model = jline.lang.Network.cyclicPs(N, D);
            elseif nargin == 3
                S = JLINE.from_line_matrix(varargin{3});
                model = jline.lang.Network.cyclicPs(N, D, S);
            else
                error('Invalid number of arguments for cyclicPs.');
            end
            model = JNetwork(model);
        end

        function model = cyclicFcfs(varargin)
            N = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            if nargin == 2
                model = jline.lang.Network.cyclicFcfs(N, D);
            elseif nargin == 3
                S = JLINE.from_line_matrix(varargin{3});
                model = jline.lang.Network.cyclicFcfs(N, D, S);
            else
                error('Invalid number of arguments for cyclicFcfs.');
            end
            model = JNetwork(model);
        end

        function model = cyclicPsInf(varargin)
            N = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            Z = JLINE.from_line_matrix(varargin{3});
            if nargin == 3
                model = jline.lang.Network.cyclicPsInf(N, D, Z);
            elseif nargin == 4
                S = JLINE.from_line_matrix(varargin{4});
                model = jline.lang.Network.cyclicPsInf(N, D, Z, S);
            else
                error('Invalid number of arguments for cyclicPsInf.');
            end
            model = JNetwork(model);
        end

        function model = cyclicFcfsInf(varargin)
            N = JLINE.from_line_matrix(varargin{1});
            D = JLINE.from_line_matrix(varargin{2});
            Z = JLINE.from_line_matrix(varargin{3});
            if nargin == 3
                model = jline.lang.Network.cyclicFcfsInf(N, D, Z);
            elseif nargin == 4
                S = JLINE.from_line_matrix(varargin{4});
                model = jline.lang.Network.cyclicFcfsInf(N, D, Z, S);
            else
                error('Invalid number of arguments for cyclicFcfsInf.');
            end
            model = JNetwork(model);
        end

        function model = cyclic(N, D, strategy, S)
            N = JLINE.from_line_matrix(N);
            D = JLINE.from_line_matrix(D);
            S = JLINE.from_line_matrix(S);
            model = jline.lang.Network.cyclic(N, D, strategy, S);
            model = JNetwork(model);
        end

        function P = serialRouting(varargin)
            if nargin == 1 && iscell(varargin{1})
                % Case: serialRouting(List<Node>)
                nodes = varargin{1};
                P = jline.lang.Network.serialRouting(nodes);

            elseif nargin == 1 && isa(varargin{1}, 'Node')
                % Case: serialRouting(Node...)
                nodeObjs = cellfun(@(n) n.obj, varargin, 'UniformOutput', false);
                P = jline.lang.Network.serialRouting(nodeObjs);

            elseif nargin == 2
                a = varargin{1};
                b = varargin{2};
                if isa(a, 'JobClass') && iscell(b)
                    % Case: serialRouting(JobClass, List<Node>)
                    P = jline.lang.Network.serialRouting(a.obj, b);

                elseif isa(a, 'JobClass') && all(cellfun(@(x) isa(x,'Node'), b))
                    % Case: serialRouting(JobClass, Node...)
                    nodeObjs = cellfun(@(n) n.obj, b, 'UniformOutput', false);
                    P = jline.lang.Network.serialRouting(a.obj, nodeObjs);

                elseif iscell(a) && iscell(b)
                    % Case: serialRouting(List<JobClass>, List<Node>)
                    P = jline.lang.Network.serialRouting(a, b);

                else
                    error('Unsupported serialRouting input pattern.');
                end

            elseif nargin >= 2 && isa(varargin{1}, 'JobClass')
                jobClass = varargin{1};
                nodes = varargin(2:end);
                nodeObjs = cellfun(@(n) n.obj, nodes, 'UniformOutput', false);
                P = jline.lang.Network.serialRouting(jobClass.obj, nodeObjs);

            elseif nargin >= 2 && all(cellfun(@(x) isa(x, 'Node'), varargin))
                % Generic case: serialRouting(Node1, Node2, ..., NodeN)
                nodeObjs = cellfun(@(n) n.obj, varargin, 'UniformOutput', false);
                P = jline.lang.Network.serialRouting([nodeObjs{:}]);

            elseif nargin == 2 && iscell(varargin{1}) && iscell(varargin{2})
                % Possibly a fallback for serialRouting(jobClasses, nodes)
                P = jline.lang.Network.serialRouting(varargin{1}, varargin{2});

            else
                error('Invalid arguments to serialRouting.');
            end
        end


    end
end
