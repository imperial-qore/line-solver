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
        admissionProb;  % q-LRU admission probability on a miss (1.0 = always admit)
        popularity;
        nLevels;
        itemLevelCap;
        items;
        accessProb;
        graph;
        totalCacheCapacity;           % sum(itemLevelCap)
        retrievalSystemCapacity;      % 0 until setRetrievalSystem is called; nitems-totalCacheCapacity otherwise
        retrievalSystemQueueIndices;  % containers.Map jobinClassIdx -> [queue node indices]
        retrievalClassIndices;        % set of indices of retrieval classes (used in afterEventCache READ)
        retrievalRoutingEntries;      % cell array of [fromCls,toCls,srcNode,dstNode,prob] routing tuples;
                                      % link() injects these into the routing matrix P (the auto-generated
                                      % retrieval classes are not part of the user-supplied P).
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
                self.totalCacheCapacity = sum(itemLevelCap);
                self.retrievalSystemCapacity = 0;  % no retrieval system by default
                if self.totalCacheCapacity > nitems
                    line_error(mfilename,sprintf('The number of items is smaller than the capacity of %s.',name));
                end
                self.retrievalSystemQueueIndices = containers.Map('KeyType','int32','ValueType','any');
                self.retrievalClassIndices = [];
                self.retrievalRoutingEntries = {};
                self.replacestrategy = replStrat;
                self.admissionProb = 1.0; % default: always admit on a miss (overridden for q-LRU)
                %probHit = min(sum(itemLevelCap)/nitems,1.0); % initial estimate of hit probability
                %self.setResultHitProb(probHit);
                %self.setResultMissProb(1-probHit);
                self.server =  CacheClassSwitcher(classes, nitems, itemLevelCap); % replace Server created by Queue
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
            self.server.actualDelayedHitProb = sparse([]);
            self.server.actualHitProbList = sparse([]);
            self.server.actualItemProb = sparse([]);
            self.server.actualResidT = sparse([]);
        end

        function self = setResultResidT(self, actualResidT)
            self.server.actualResidT = actualResidT;
        end

        function p = getResidT(self)
            p = full(self.server.actualResidT);
        end

        function tc = getTotalCacheCapacity(self)
            tc = self.totalCacheCapacity;
        end

        function rc = getRetrievalSystemCapacity(self)
            rc = self.retrievalSystemCapacity;
        end

        function m = getRetrievalClasses(self)
            m = self.server.retrievalClasses;
        end

        function idx = getRetrievalClassIndices(self)
            idx = self.retrievalClassIndices;
        end

        function q = getRetrievalSystemQueueIndicesFor(self, jobinClassIdx)
            % q = getRetrievalSystemQueueIndicesFor(jobinClassIdx)
            % Return node indices of the queues comprising the retrieval system for the
            % given (0-indexed) arrival class, or [] if no retrieval system is set.
            if isKey(self.retrievalSystemQueueIndices, int32(jobinClassIdx))
                q = self.retrievalSystemQueueIndices(int32(jobinClassIdx));
            else
                q = [];
            end
        end

        function setRetrievalClass(self, jobinClass, joboutClass, item)
            % SETRETRIEVALCLASS(jobinClass, joboutClass, item)
            % item is 1-based (MATLAB convention).
            self.server.retrievalClasses(item, jobinClass.index) = joboutClass.index;
        end
                
        function self = setResultHitProb(self, actualHitProb)
            self.server.actualHitProb = actualHitProb;
        end
        
        function self = setResultMissProb(self, actualMissProb)
            self.server.actualMissProb = actualMissProb;
        end

        function self = setResultDelayedHitProb(self, actualDelayedHitProb)
            % SETRESULTDELAYEDHITPROB  Per-class delayed-hit fraction
            % (requests arriving for an item whose fetch is already in
            % progress in the retrieval system). Zero for caches without a
            % retrieval system.
            self.server.actualDelayedHitProb = actualDelayedHitProb;
        end

        function p = getHitRatio(self)
            % GETHITRATIO  Actual (true) hit fraction per class: the item is
            % resident in the cache. Delayed hits are reported separately by
            % getDelayedHitRatio.
            p = full(self.server.actualHitProb);
        end

        function p = getMissRatio(self)
            p = full(self.server.actualMissProb);
        end

        function p = getDelayedHitRatio(self)
            % GETDELAYEDHITRATIO  Actual delayed-hit fraction per class
            % (empty/zero when the cache has no retrieval system).
            if isprop(self.server, 'actualDelayedHitProb')
                p = full(self.server.actualDelayedHitProb);
            else
                p = [];
            end
        end

        function self = setResultHitProbList(self, actualHitProbList)
            % SETRESULTHITPROBLIST  Per-class, per-list (per-level) hit
            % fraction matrix [classes x lists]; rows sum to getHitRatio.
            self.server.actualHitProbList = actualHitProbList;
        end

        function p = getHitRatioByList(self)
            % GETHITRATIOBYLIST  Per-class, per-list hit fraction matrix
            % [classes x lists]; empty when not computed by the solver.
            if isprop(self.server, 'actualHitProbList')
                p = full(self.server.actualHitProbList);
            else
                p = [];
            end
        end

        function self = setResultItemProb(self, actualItemProb)
            % SETRESULTITEMPROB  Per-item occupancy matrix [items x (lists+1)];
            % column 1 = miss (item not cached), columns 2..end = per-list.
            self.server.actualItemProb = actualItemProb;
        end

        function p = getItemProb(self)
            % GETITEMPROB  Per-item occupancy matrix [items x (lists+1)]: column 1
            % is the miss probability, columns 2..end the per-list probabilities;
            % empty when not computed by the solver.
            if isprop(self.server, 'actualItemProb')
                p = full(self.server.actualItemProb);
            else
                p = [];
            end
        end
        
        function setHitClass(self, jobinclass, joboutclass)
            % SETHITCLASS(JOBINCLASS, JOBOUTCLASS)
            
            self.server.hitClass(jobinclass.index) = joboutclass.index;
        end
        
        function setMissClass(self, jobinclass, joboutclass)
            % SETMISSCLASS(JOBINCLASS, JOBOUTCLASS)
            
            self.server.missClass(jobinclass.index) = joboutclass.index;
        end
        
       
        function self = removeJobClass(self, jobclass)
            % SELF = REMOVEJOBCLASS(JOBCLASS)
            %
            % Reject class removal: the cache item state is indexed by class
            % and lives in the model state, not only in this node, so it
            % cannot be re-indexed here. Matches Cache.removeJobClass in the
            % JAR and Cache.remove_job_class in python.

            line_error(mfilename,'Cannot dynamically remove classes in models with caches. You need to re-instantiate the model.');
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

        function setAdmissionProb(self, q)
            % SETADMISSIONPROB(Q)
            % Probability q in [0,1] of admitting a missed item into the cache
            % (q-LRU). Only used when the replacement strategy is QLRU.
            if q < 0 || q > 1
                line_error(mfilename,'The admission probability q must lie in [0,1].');
            end
            self.admissionProb = q;
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

        function addRetrievalRoutingEntry(self, fromCls, toCls, srcNode, dstNode, prob, allowZero)
            % ADDRETRIEVALROUTINGENTRY(fromCls, toCls, srcNode, dstNode, prob, allowZero)
            % Register a routing edge for an auto-generated retrieval class. link()
            % injects all such entries into the routing matrix P. Entries are 1-based
            % class indices and 1-based node indices; prob is the routing probability.
            % Later entries override earlier ones for the same (fromCls,toCls,src,dst).
            % With allowZero=true an explicit prob==0 entry is recorded so it can
            % override (delete) a default edge inherited from the read class; the
            % internal broadcast path keeps allowZero=false and drops zero edges.
            if nargin < 7
                allowZero = false;
            end
            if prob < 0 || (prob == 0 && ~allowZero)
                return
            end
            self.retrievalRoutingEntries{end+1} = [fromCls, toCls, srcNode, dstNode, prob];
        end

        function setItemRoutingProbability(self, jobinClass, item, source, dest, probability)
            % SETITEMROUTINGPROBABILITY(jobinClass, item, source, dest, probability)
            % Probability of routing the retrieval class for `item` between two nodes of
            % the retrieval system. `source`/`dest` are either a retrieval queue or the
            % cache itself: pass the cache as `source` for a cache->queue entry, or as
            % `dest` for a queue->cache exit.
            rClassIdx = self.server.retrievalClasses(item, jobinClass.index);
            if rClassIdx <= 0
                line_error(mfilename,'No retrieval class defined for the given class/item; call setRetrievalSystem first.');
            end
            self.addRetrievalRoutingEntry(rClassIdx, rClassIdx, source.index, dest.index, probability, true);
        end

        function setItemRoutingProb(self, jobinClass, item, source, dest, probability)
            % Short alias for setItemRoutingProbability.
            self.setItemRoutingProbability(jobinClass, item, source, dest, probability);
        end

        function setRetrievalSystem(self, jobinClass, missClass, queues)
            % SETRETRIEVALSYSTEM(jobinClass, missClass, queues)
            %
            % Initialise the retrieval system through which a request that misses the cache
            % is fetched. The request switches to a per-item retrieval class that circulates
            % the queues and returns to the cache, where the returning READ logs it as a
            % miss (switching into `missClass`). getResidT reports the per-class
            % queueing time in the retrieval sub-network.
            %
            % Arguments:
            %   jobinClass      arrival JobClass that can route through the retrieval system
            %   missClass       JobClass into which a completed retrieval transitions
            %   queues          single Queue or array/cell of Queue nodes comprising the system
            %
            % Routing and service are NOT passed here; they are taken from the read class:
            %   - service: the read class's service distribution at each queue. Call
            %     queue.setService(jobinClass, ...) beforehand; override per item with
            %     queue.setItemServiceRate(cache, jobinClass, item, rate).
            %   - routing: the read class's routing among the retrieval queues drawn in the
            %     top-level routing matrix P; override per item with
            %     setItemQueueEntryProbability / setItemRoutingProbability / setItemQueueExitProbability.

            % --- normalise inputs ---
            if isa(queues, 'Queue')
                queueArr = {queues};
            elseif iscell(queues)
                queueArr = queues;
            elseif isnumeric(queues)
                line_error(mfilename,'queues must be a Queue or a cell/array of Queues.');
            else
                queueArr = num2cell(queues);
            end
            nQueues = numel(queueArr);
            nItems = self.items.nitems;
            self.retrievalSystemCapacity = nItems - self.totalCacheCapacity;

            if nQueues == 0
                line_error(mfilename,'Retrieval system cannot be initialised with no stations.');
            end

            % inherit the read class's service distribution at each queue as the per-item default
            serviceDistByQueue = cell(1, nQueues);
            for q = 1:nQueues
                sp = queueArr{q}.serviceProcess;
                if numel(sp) >= jobinClass.index && ~isempty(sp{jobinClass.index})
                    serviceDistByQueue{q} = sp{jobinClass.index};
                else
                    line_error(mfilename, sprintf(['No service distribution for the read class at queue "%s"; ' ...
                        'call queue.setService(readClass, ...) before setRetrievalSystem.'], queueArr{q}.name));
                end
            end

            % --- record queue node indices ---
            queueIdxs = zeros(1, nQueues);
            for q = 1:nQueues
                queueIdxs(q) = queueArr{q}.index;
            end
            self.retrievalSystemQueueIndices(int32(jobinClass.index-1)) = queueIdxs;

            % --- create one retrieval class per item ---
            retrievalList = cell(1, nItems);
            for i = 1:nItems
                if isa(jobinClass, 'ClosedClass')
                    refStation = jobinClass.refstat;
                    retrievalList{i} = ClosedClass(self.model, [jobinClass.name '_retrievalClass_' num2str(i)], 0, refStation, 0);
                else
                    retrievalList{i} = OpenClass(self.model, [jobinClass.name '_retrievalClass_' num2str(i)], 0);
                end
                self.retrievalClassIndices(end+1) = retrievalList{i}.index;
            end

            % --- per-item retrieval class setup ---
            % Each item's retrieval class reads item i (one-hot popularity) and, on the
            % returning READ, is logged as a miss. Its routing through the queues is NOT set
            % here: it is inherited at link() from the read class's routing in P, overridable
            % per item via the setItem* methods.
            for i = 1:nItems
                rClass = retrievalList{i};

                % switch arrival jobinClass -> retrieval class for item i
                self.setRetrievalClass(jobinClass, rClass, i);

                % the retrieval class triggers a read of item i (one-hot popularity)
                itemPopularity = zeros(1, nItems);
                itemPopularity(i) = 1.0;
                self.popularity{self.items.index, rClass.index} = DiscreteSampler(itemPopularity);

                % on the returning READ the retrieval is logged as a miss
                self.setMissClass(rClass, missClass);

                for sourceQueueIdx = 1:nQueues
                    sourceQueue = queueArr{sourceQueueIdx};

                    % service for the retrieval class at this queue (inherited read-class service)
                    sourceQueue.setService(rClass, serviceDistByQueue{sourceQueueIdx}.copy());

                    % at most one retrieval in flight at a time for this class
                    if length(sourceQueue.classCap) < rClass.index
                        sourceQueue.classCap((length(sourceQueue.classCap)+1):rClass.index) = Inf;
                    end
                    sourceQueue.classCap(rClass.index) = 1;
                end
            end
        end
    end
end
