classdef LayeredNetwork < Model & Ensemble
    % A layered queueing network model.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties (Hidden)
        usedFeatures; % cell with structures of booleans listing the used classes
        % it must be accessed via getUsedLangFeatures
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
            % SELF = LAYEREDNETWORK(NAME, FILENAME)

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
        name = getNodeName(self,node)
        node = getNodeByName(self,name)
        [names,hostnames,tasknames,entrynames,actnames] = getNodeNames(self)
        proc = getNodeHost(self,node)
        task = getNodeTask(self,node)
        type = getNodeType(self,node)

        writeSRVN(self,filename);
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
        function self = setUsedFeatures(self,e,className)
            % SELF = SETUSEDFEATURES(SELF,E,CLASSNAME)

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

        function [lambda,D,N,Z,mu,S]= getProductFormParameters(self)
            % [LAMBDA,D,N,Z,MU,S]= GETPRODUCTFORMPARAMETERS()

            ensemble = self.getEnsemble;
            for e=1:length(ensemble)
                [lambda{e,1},D{e,1},N{e,1},Z{e,1},mu{e,1},S{e,1}] = ensemble{e}.getProductFormParameters;
            end
        end

        function [lambda,D,N,Z,mu,S]= getProductFormChainParameters(self)
            % [LAMBDA,D,N,Z,MU,S]= GETPRODUCTFORMCHAINPARAMETERS()

            ensemble = self.getEnsemble;
            for e=1:length(ensemble)
                [lambda{e,1},D{e,1},N{e,1},Z{e,1},mu{e,1},S{e,1}] = ensemble{e}.getProductFormChainParameters;
            end
        end

        function nnodes = getNumberOfStatefulNodes(self)
            if isempty(self.ensemble)
                self.getEnsemble();
            end
            nnodes = 0; % delay and queue
            for e=1:length(self.ensemble)
                nnodes = nnodes + self.ensemble{e}.getNumberOfStatefulNodes();
            end
        end

    end

    methods (Static)
        function myLN = readXML(filename, verbose)
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
