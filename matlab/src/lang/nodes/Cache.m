classdef Cache < StatefulNode
    % Multi-level cache node with hit/miss class switching
    %
    % Models cache memory systems with multiple levels and replacement strategies.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        cap;
        schedPolicy;
        schedStrategy;
        replacestrategy;
        popularity;
        nLevels;
        itemLevelCap;
        items;
        accessProb;
        graph;
    end
    
    methods
        %Constructor
        function self = Cache(model, name, nitems, itemLevelCap, replStrat, graph)
            % CACHE Create a Cache node instance
            %
            % @brief Creates a Cache node with configurable levels and replacement strategy
            % @param model Network model to add the cache to
            % @param name String identifier for the cache node
            % @param nitems Total number of cacheable items
            % @param itemLevelCap Vector specifying capacity of each cache level
            % @param replStrat Replacement strategy (LRU, FIFO, Random, etc.)
            % @param graph Optional graph structure for cache hierarchy
            % @return self Cache instance configured for the given model
            %
            % The constructor creates a multi-level cache with the specified
            % total items and per-level capacities. The replacement strategy
            % determines how items are evicted when cache levels become full.
            
            self@StatefulNode(name);
            if model.isMatlabNative()
                if ~exist('itemLevelCap','var')
                    levels = 1;
                end
                classes = model.getClasses();
                self.input = Buffer(classes);
                self.output = Dispatcher(classes);
                self.schedPolicy = SchedStrategyType.NP;
                self.schedStrategy = SchedStrategy.FCFS;
                self.items = ItemSet(model, [name,'_','Items'], nitems, self);
                self.nLevels = nnz(itemLevelCap);
                self.cap = Inf; % job capacity
                self.accessProb = {};
                self.itemLevelCap = itemLevelCap; % item capacity
                if sum(itemLevelCap) > nitems
                    line_error(mfilename,sprintf('The number of items is smaller than the capacity of %s.',name));
                end
                self.replacestrategy = replStrat;
                %probHit = min(sum(itemLevelCap)/nitems,1.0); % initial estimate of hit probability
                %self.setResultHitProb(probHit);
                %self.setResultMissProb(1-probHit);
                self.server =  CacheClassSwitcher(classes, self.nLevels, itemLevelCap); % replace Server created by Queue
                self.popularity = {};
                self.setModel(model);
                self.model.addNode(self);
                if nargin<6
                    self.graph = [];
                else
                    self.graph = graph;
                end
            elseif model.isJavaNative()
                self.setModel(model);
                if nargin<6 || isempty(graph)
                    self.obj = jline.lang.nodes.Cache(model.obj, name, nitems, itemLevelCap, replStrat);
                else
                    self.obj = jline.lang.nodes.Cache(model.obj, name, nitems, itemLevelCap, replStrat, graph);
                end
                self.index = model.obj.getNodeIndex(self.obj);
            end            
        end
        
        %         function setMissTime(self, distribution)
        % SETMISSTIME(DISTRIBUTION)
        
        %             itemclass = self.items;
        %             self.server.serviceProcess{1, itemclass.index} = {[], ServiceStrategy.SD, distribution};
        %         end
        %
        %         function setHitTime(self, distribution, level)
        % SETHITTIME(DISTRIBUTION, LEVEL)
        
        %             itemclass = self.items;
        %             if ~exist('level','var')
        %                 levels = 2:self.nLevels;
        %             else
        %                 levels = level;
        %             end
        %             for level = levels
        %                 self.server.serviceProcess{1+level, itemclass.index} = {[], ServiceStrategy.SD, distribution};
        %             end
        %         end
        
        function self = reset(self)
            % SELF = RESET()
            %
            % Reset internal data structures when the network model is
            % reset
            self.server.actualHitProb = sparse([]);
            self.server.actualMissProb = sparse([]);
        end
                
        function self = setResultHitProb(self, actualHitProb)
            self.server.actualHitProb = actualHitProb;
        end
        
        function self = setResultMissProb(self, actualMissProb)
            self.server.actualMissProb = actualMissProb;
        end
        
        function p = getHitRatio(self)
            p = full(self.server.actualHitProb);
        end

        function p = getMissRatio(self)
            p = full(self.server.actualMissProb);
        end
        
        function setHitClass(self, jobinclass, joboutclass)
            % SETHITCLASS(JOBINCLASS, JOBOUTCLASS)
            
            self.server.hitClass(jobinclass.index) = joboutclass.index;
        end
        
        function setMissClass(self, jobinclass, joboutclass)
            % SETMISSCLASS(JOBINCLASS, JOBOUTCLASS)
            
            self.server.missClass(jobinclass.index) = joboutclass.index;
        end
        
       
        function setRead(self, jobclass, distribution)
            % SETREAD(JOBCLASS, DISTRIBUTION)
            
            itemclass = self.items;
            if distribution.isDiscrete                
                self.popularity{itemclass.index, jobclass.index} = distribution.copy;
                if self.popularity{itemclass.index, jobclass.index}.support(2) ~= itemclass.nitems
                    line_error(mfilename,sprintf('The reference model is defined on a number of items different from the ones used to instantiate %s.',self.name));
                end
                switch class(distribution)
                    case 'Zipf'
                        self.popularity{itemclass.index, jobclass.index}.setParam(2, 'n', itemclass.nitems);
                end
                %                self.probselect(itemclass.index, jobclass.index) = probselect;
            else
                line_error(mfilename,'A discrete popularity distribution is required.');
            end
        end
        
        function setReadItemEntry(self, jobclass, popularity, cardinality)
            % SETREAD(JOBCLASS, DISTRIBUTION)
            
            if popularity.isDiscrete
                
                self.popularity{jobclass.index} = popularity.copy;
                switch class(popularity)
                    case 'Zipf'
                        self.popularity{jobclass.index}.setParam(2, 'n', cardinality);
                end
                
            else
                line_error(mfilename,'A discrete popularity distribution is required.');
            end
        end
        function setAccessProb(self, R)
            % SETACCESSCOSTS(R)
            
            self.accessProb = R;
        end
        
        
        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end
                
        function hitClass = getHitClass(self)
            % HITCLASS = GETHITCLASS
            %
            % For an incoming job of class r, HITCLASS(r) is the new class
            % of that job after a hit
            
            hitClass = self.server.hitClass;
        end
        
        function missClass = getMissClass(self)
            % MISSCLASS = GETMISSCLASS
            %
            % For an incoming job of class r, MISSCLASS(r) is the new class
            % of that job after a miss
            
            missClass = self.server.missClass;
        end
    end
end
