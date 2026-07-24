classdef LayeredNetwork < Model & Ensemble
    % LayeredNetwork Hierarchical software performance modeling framework
    %
    % LayeredNetwork implements Layered Queueing Networks (LQN) for modeling
    % hierarchical software systems with clients, application servers, and
    % resource layers. It supports modeling of complex software architectures
    % with tasks, entries, activities, and their interactions across multiple
    % system layers.
    %
    % @brief Layered queueing network for hierarchical software performance models
    %
    % Key characteristics:
    % - Hierarchical multi-layer architecture modeling
    % - Software-centric performance analysis
    % - Task and entry abstraction levels
    % - Activity-based detailed modeling
    % - Host resource modeling
    % - Client-server interaction patterns
    %
    % LQN model components:
    % - Hosts: Physical or logical processing resources
    % - Tasks: Software processes or services
    % - Entries: Service request entry points
    % - Activities: Detailed task execution steps
    % - Reference tasks: External workload generators
    %
    % LayeredNetwork is used for:
    % - Software performance engineering
    % - Multi-tier application modeling
    % - Microservice architecture analysis
    % - Distributed system performance evaluation
    % - Capacity planning for software systems
    %
    % Example:
    % @code
    % lqn = LayeredNetwork('WebApp');
    % client = Host(lqn, 'ClientTier', Inf);
    % server = Host(lqn, 'ServerTier', 4);
    % web_task = Task(lqn, 'WebServer', 1, server);
    % request_entry = Entry(lqn, 'ProcessRequest', web_task);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Hidden)
        lsn;
        usedFeatures; % cell with structures of booleans listing the used classes
        % it must be accessed via getUsedLangFeatures
        obj = []; % jline.lang.layered.LayeredNetwork mirror (set by Java-backed solvers)
    end

    properties
        hosts = [];
        tasks = [];
        reftasks = [];
        activities = [];
        entries = [];
    end

    methods
        %public methods, including constructor

        function self = reset(self, isHard)
            if nargin<2
                isHard = false;
            end
            self.ensemble = {};
            if isHard
                self.hosts = {};
                self.tasks = {};
                self.reftasks = {};
                self.entries = {};
                self.activities = {};
            end
        end

        % constructor
        function self = LayeredNetwork(name, filename)
            % LAYEREDNETWORK Create a layered queueing network model
            %
            % @brief Creates a LayeredNetwork instance for hierarchical modeling
            % @param name String identifier for the layered model
            % @param filename Optional filename for model import/export
            % @return self LayeredNetwork instance ready for hierarchical modeling

            self@Ensemble({})
            if nargin<1 %~exist('name','var')
                [~,name]=fileparts(lineTempName);
            end
            name = char(name);
            self@Model(name);
            self.ensemble = {};
            self.hosts = {};
            self.tasks = {};
            self.reftasks = {};
            self.entries = {};
            self.activities = {};

            if nargin>=2 %exist('filename','var')
                self = LayeredNetwork.parseXML(filename, false);
            end
        end


        function sn = summary(self)
            % sn = SUMMARY()

            sn = self.getStruct;
        end

        plot(self, showTaskGraph)
        plotGraph(self, useNodes)
        plotGraphSimple(self, useNodes)
        plotTaskGraph(self, useNodes)
    end

    methods
        idx = getNodeIndex(self,node)
        node = getNodeByName(self,name)
        [names,hostnames,tasknames,entrynames,actnames] = getNodeNames(self)        

        writeXML(self,filename,useAbstractNames);
    end

    methods

        LQN = getStruct(self);

        function E = getNumberOfLayers(self)
            % E = GETNUMBEROFLAYERS()

            E = getNumberOfModels(self);
        end

        function E = getNumberOfModels(self)
            % E = GETNUMBEROFMODELS()

            if isempty(self.ensemble)
                self.ensemble = getEnsemble(self);
            end
            E = length(self.ensemble);
        end

        function layers = getLayers(self)
            % LAYERS = GETLAYERS()

            layers = getEnsemble(self);
        end

        % setUsedFeatures : records that a certain language feature has been used
        function self = setUsedLangFeature(self,e,className)
            % SELF = SETUSEDLANGFEATURE(SELF,E,CLASSNAME)

            self.usedFeatures{e}.setTrue(className);
        end

        function self = initUsedFeatures(self)
            % SELF = INITUSEDFEATURES()

            for e=1:getNumberOfModels(self)
                self.usedFeatures{e} = SolverFeatureSet;
            end
        end

        function usedFeatures = getUsedLangFeatures(self)
            % USEDFEATURES = GETUSEDLANGFEATURES()

            E = getNumberOfLayers(self);
            usedFeatures = cell(1,E);
            for e=1:E
                usedFeatures{e} = self.ensemble{e}.getUsedLangFeatures;
            end
            self.usedFeatures = usedFeatures;
        end

        function view(self)
            jlqnmodel = JLINE.from_line_layered_network(self);
            jlqnmodel.view();
        end

        function result = nodeIndex(self, varargin)
            % NODEINDEX Kotlin-style alias for getNodeIndex
            result = self.getNodeIndex(varargin{:});
        end

        function result = nodeByName(self, varargin)
            % NODEBYNAME Kotlin-style alias for getNodeByName
            result = self.getNodeByName(varargin{:});
        end

        function result = nodeNames(self, varargin)
            % NODENAMES Kotlin-style alias for getNodeNames
            result = self.getNodeNames(varargin{:});
        end

        function result = struct(self, varargin)
            % STRUCT Kotlin-style alias for getStruct
            result = self.getStruct(varargin{:});
        end

        function result = numberOfLayers(self, varargin)
            % NUMBEROFLAYERS Kotlin-style alias for getNumberOfLayers
            result = self.getNumberOfLayers(varargin{:});
        end

        function result = numberOfModels(self, varargin)
            % NUMBEROFMODELS Kotlin-style alias for getNumberOfModels
            result = self.getNumberOfModels(varargin{:});
        end

        function result = layers(self, varargin)
            % LAYERS Kotlin-style alias for getLayers
            result = self.getLayers(varargin{:});
        end

        function result = usedLangFeatures(self, varargin)
            % USEDLANGFEATURES Kotlin-style alias for getUsedLangFeatures
            result = self.getUsedLangFeatures(varargin{:});
        end

        function sanitize(self)
            % SANITIZE()
            % Validates the LayeredNetwork configuration.
            % Ensures that if entries are defined, activities are also defined to serve those entries.

            numEntries = length(self.entries);
            numActivities = length(self.activities);

            if numEntries > 0 && numActivities == 0
                msg = sprintf('LayeredNetwork ''%s'' has %d entry(ies) but no activities. Entries must be bound to activities to form a valid LQN model. Use activity.boundTo(entry) to establish the binding.', self.name, numEntries);
                line_error(mfilename, msg);
            end
        end

        % Aggregate flat interface (Network-compatible) used by SolverENV.
        % An LQN is exposed as the block-diagonal union of its layer
        % networks: aggregate station/class/node counts are the sums over
        % layers, and matrices are laid out block-diagonally in the exact
        % ordering produced by SolverLN.getTranAvg. layerBlocks returns the
        % per-layer cumulative offsets and sizes that all these methods and
        % the SolverENV labeling share as a single source of truth.
        function [Roff, Coff, Msz, Ksz] = layerBlocks(self)
            % [ROFF,COFF,MSZ,KSZ] = LAYERBLOCKS()
            % Roff(e)/Coff(e): row/col offset (0-based) of layer e's block;
            % Msz(e)/Ksz(e): number of stations/classes of layer e.
            if isempty(self.ensemble)
                self.ensemble = getEnsemble(self);
            end
            E = length(self.ensemble);
            Msz = zeros(1,E); Ksz = zeros(1,E);
            for e=1:E
                Msz(e) = self.ensemble{e}.getNumberOfStations;
                Ksz(e) = self.ensemble{e}.getNumberOfClasses;
            end
            Roff = [0, cumsum(Msz(1:end-1))];
            Coff = [0, cumsum(Ksz(1:end-1))];
        end

        function M = getNumberOfStations(self)
            % M = GETNUMBEROFSTATIONS()
            % Aggregate (sum over layers) station count.
            [~,~,Msz] = self.layerBlocks;
            M = sum(Msz);
        end

        function K = getNumberOfClasses(self)
            % K = GETNUMBEROFCLASSES()
            % Aggregate (sum over layers) class count.
            [~,~,~,Ksz] = self.layerBlocks;
            K = sum(Ksz);
        end

        function N = getNumberOfNodes(self)
            % N = GETNUMBEROFNODES()
            % Aggregate (sum over layers) node count.
            if isempty(self.ensemble)
                self.ensemble = getEnsemble(self);
            end
            N = 0;
            for e=1:length(self.ensemble)
                N = N + self.ensemble{e}.getNumberOfNodes;
            end
        end

        function N = getNumberOfStatefulNodes(self)
            % N = GETNUMBEROFSTATEFULNODES()
            % Aggregate (sum over layers) stateful-node count. Used by
            % Environment.addStage to validate stage compatibility.
            if isempty(self.ensemble)
                self.ensemble = getEnsemble(self);
            end
            N = 0;
            for e=1:length(self.ensemble)
                N = N + self.ensemble{e}.getNumberOfStatefulNodes;
            end
        end

        function [Qt, Ut, Tt] = getTranHandles(self)
            % [QT,UT,TT] = GETTRANHANDLES()
            % Block-diagonal aggregate transient handles over the layer
            % networks. Off-diagonal cells are left empty (no cross-layer
            % station/class pairing); SolverLN.getTranAvg reproduces the same
            % block-diagonal layout so the handles only fix the M x K shape.
            [Roff, Coff, Msz, Ksz] = self.layerBlocks;
            M = sum(Msz); K = sum(Ksz);
            Qt = cell(M,K); Ut = cell(M,K); Tt = cell(M,K);
            for e=1:length(self.ensemble)
                [Qe, Ue, Te] = self.ensemble{e}.getTranHandles;
                for i=1:Msz(e)
                    for r=1:Ksz(e)
                        Qt{Roff(e)+i, Coff(e)+r} = Qe{i,r};
                        Ut{Roff(e)+i, Coff(e)+r} = Ue{i,r};
                        Tt{Roff(e)+i, Coff(e)+r} = Te{i,r};
                    end
                end
            end
        end

        function initFromMarginal(self, n, options)
            % INITFROMMARGINAL(N, OPTIONS)
            % Split the aggregate (M x K) marginal queue-length matrix N into
            % per-layer blocks and warm-start each layer network. This is the
            % state-continuity mechanism used by SolverENV across environment
            % switches.
            if nargin<3
                options = Solver.defaultOptions;
            end
            [Roff, Coff, Msz, Ksz] = self.layerBlocks;
            for e=1:length(self.ensemble)
                block = n(Roff(e)+1:Roff(e)+Msz(e), Coff(e)+1:Coff(e)+Ksz(e));
                self.ensemble{e}.initFromMarginal(block, options);
            end
        end

        % Getter methods for API consistency with Java/Python
        function val = getHosts(obj)
            % GETHOSTS Get the hosts/processors in the network
            val = obj.hosts;
        end

        function val = getTasks(obj)
            % GETTASKS Get the tasks in the network
            val = obj.tasks;
        end

        function val = getEntries(obj)
            % GETENTRIES Get the entries in the network
            val = obj.entries;
        end

        function val = getActivities(obj)
            % GETACTIVITIES Get the activities in the network
            val = obj.activities;
        end

    end

    methods (Static)
        function myLN = readXML(filename, verbose)
            if nargin < 2
                verbose = false;
            end
            myLN = LayeredNetwork.parseXML(filename, verbose);
        end
        
        function myLN = load(filename, verbose)
            if nargin < 2
                verbose = false;
            end
            myLN = LayeredNetwork.parseXML(filename, verbose);
        end
        
        myLN = parseXML(filename, verbose)

        function myLN = fromNetwork(model)
            myLN = QN2LQN(model);
        end
    end
end
