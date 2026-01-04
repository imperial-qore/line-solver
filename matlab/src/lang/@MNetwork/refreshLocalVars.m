function nvars = refreshLocalVars(self)
% NVARS = REFRESHLOCALVARS()

R = self.getNumberOfClasses;
nvars = zeros(self.getNumberOfNodes, 2*R+1);
nodeparam = cell(self.getNumberOfNodes, 1);
rtnodes = self.sn.rtnodes;
% Draft SPN code:
% isp = [];
% ist = [];
% nodeToPlace = zeros(1, self.getNumberOfNodes);
% nodeToTransition = zeros(1, self.getNumberOfNodes);
% for ind=1:self.getNumberOfNodes
%     node = self.getNodeByIndex(ind);
%     switch class(node)
%     case 'Place'
%         isp = [isp, ind];
%         nodeToPlace(ind) = length(isp);
%     case 'Transition'
%         ist = [ist, ind];
%         nodeToTransition(ind) = length(ist);
%     end
% end
sn = self.sn;

for ind=1:self.getNumberOfNodes
    node = self.getNodeByIndex(ind);
    switch class(node)
        case 'Cache'
            nvars(ind,2*R+1) = sum(node.itemLevelCap);
            nodeparam{ind} = struct();
            nodeparam{ind}.nitems = 0;
            nodeparam{ind}.accost = node.accessProb;
            for r=1:self.getNumberOfClasses
                if length(node.popularity) >= r && isa(node.popularity{r}, 'Distribution') && ~node.popularity{r}.isDisabled
                    nodeparam{ind}.nitems = max(nodeparam{ind}.nitems,node.popularity{r}.support(2));
                end
            end
            nodeparam{ind}.itemcap = node.itemLevelCap;
            nodeparam{ind}.pread = cell(1,self.getNumberOfClasses);
            for r=1:self.getNumberOfClasses
                if length(node.popularity) < r || ~isa(node.popularity{r}, 'Distribution') || node.popularity{r}.isDisabled
                    nodeparam{ind}.pread{r} = NaN;
                else
                    nodeparam{ind}.pread{r} = node.popularity{r}.evalPMF(1:nodeparam{ind}.nitems);
                end
            end
            nodeparam{ind}.replacestrat = node.replacestrategy;
            nodeparam{ind}.hitclass = zeros(1,self.getNumberOfClasses);
            nodeparam{ind}.hitclass(1:length(node.server.hitClass)) = round(node.server.hitClass);
            nodeparam{ind}.missclass = zeros(1,self.getNumberOfClasses);
            nodeparam{ind}.missclass(1:length(node.server.missClass)) = round(node.server.missClass);
            % Store actual hit/miss probabilities from solver results (if available)
            if isprop(node.server, 'actualHitProb') && ~isempty(node.server.actualHitProb)
                nodeparam{ind}.actualhitprob = full(node.server.actualHitProb);
                nodeparam{ind}.actualmissprob = full(node.server.actualMissProb);
            end
        case 'Fork'
            nodeparam{ind}.fanOut = node.output.tasksPerLink;
        case 'Join'
            nodeparam{ind}.joinStrategy = node.input.joinStrategy;
            nodeparam{ind}.fanIn = cell(1,self.getNumberOfClasses);
            for r=1:self.getNumberOfClasses
                nodeparam{ind}.fanIn{r} = node.input.joinRequired{r};
            end
        case 'Logger'
            nodeparam{ind}.fileName = node.fileName;
            nodeparam{ind}.filePath = node.filePath;
            nodeparam{ind}.startTime = node.getStartTime;
            nodeparam{ind}.loggerName = node.getLoggerName;
            nodeparam{ind}.timestamp = node.getTimestamp;
            nodeparam{ind}.jobID = node.getJobID;
            nodeparam{ind}.jobClass = node.getJobClass;
            nodeparam{ind}.timeSameClass = node.getTimeSameClass;
            nodeparam{ind}.timeAnyClass = node.getTimeAnyClass;            
        case {'Source'}
            for r=1:self.getNumberOfClasses
                if ~iscell(node.arrivalProcess) || length(node.arrivalProcess) < r || isempty(node.arrivalProcess{r})
                    continue;
                end
                switch class(node.arrivalProcess{r})
                    case 'MAP'
                        nvars(ind,r) = nvars(ind,r) + 1;
                    case {'Replayer', 'Trace'}
                        if isempty(nodeparam{ind})
                            nodeparam{ind} = cell(1,self.getNumberOfClasses);
                        end
                        if isempty(nodeparam{ind}{r})
                            nodeparam{ind}{r} = struct();
                        end
                        nodeparam{ind}{r}.(node.arrivalProcess{r}.params{1}.paramName) = node.arrivalProcess{r}.params{1}.paramValue;
                end
            end
        case {'Queue','QueueingStation','Delay','DelayStation','Transition'}
            for r=1:self.getNumberOfClasses
                switch class(node.server.serviceProcess{r}{3})
                    case 'MAP'
                        nvars(ind,r) = nvars(ind,r) + 1;
                    case {'Replayer', 'Trace'}
                        if isempty(nodeparam{ind})
                            nodeparam{ind} = cell(1,self.getNumberOfClasses);
                        end
                        if isempty(nodeparam{ind}{r})
                            nodeparam{ind}{r} = struct();
                        end
                        nodeparam{ind}{r}.(node.server.serviceProcess{r}{3}.params{1}.paramName) = node.server.serviceProcess{r}{3}.params{1}.paramValue;
                end
                if ~isa(node,'Delay') && ~isa(node,'Transition') && ~isempty(node.setupTime) && ~isempty(node.setupTime{r})
                    if isempty(nodeparam{ind})
                        nodeparam{ind} = cell(1,self.getNumberOfClasses);
                    end
                    if isempty(nodeparam{ind}{r})
                        nodeparam{ind}{r} = struct();
                    end
                    nodeparam{ind}{r}.setupTime = node.setupTime{r}.getProcess();
                    nodeparam{ind}{r}.delayoffTime = node.delayoffTime{r}.getProcess();
                end
                if (isa(node,'Queue') || isa(node,'QueueingStation'))
                    if isempty(nodeparam{ind}) && (~isempty(node.pollingType) || ~isempty(node.switchoverTime))
                        nodeparam{ind} = cell(1,self.getNumberOfClasses);
                        for s=1:self.getNumberOfClasses
                            nodeparam{ind}{s} = struct();
                        end
                    end
                    if ~isempty(node.pollingType) && ~isempty(node.pollingType{r})
                        nodeparam{ind}{r}.pollingType = node.pollingType{r};
                        nodeparam{ind}{r}.pollingPar = node.pollingPar;
                    end
                    if ~isempty(node.switchoverTime) && ~isempty(node.switchoverTime{r})
                        if min(size(node.switchoverTime))==1
                            nodeparam{ind}{r}.switchoverTime = node.switchoverTime{r}.getProcess();
                            nodeparam{ind}{r}.switchoverProcId = ProcessType.toId(ProcessType.fromText(class(node.switchoverTime{r})));
                        else
                            for t=1:length(node.switchoverTime)
                                nodeparam{ind}{r}.switchoverTime{t} = node.switchoverTime{r,t}.getProcess();
                                nodeparam{ind}{r}.switchoverProcId(t) = ProcessType.toId(ProcessType.fromText(class(node.switchoverTime{r,t})));
                            end
                        end
                    end
                end
            end
    end

    for r=1:R
        switch sn.routing(ind,r)
            case RoutingStrategy.KCHOICES
                nodeparam{ind}{r}.k = node.output.outputStrategy{r}{3}{1};
                nodeparam{ind}{r}.withMemory = node.output.outputStrategy{r}{3}{2};
            case RoutingStrategy.WRROBIN
                nvars(ind,R+r) = nvars(ind,R+r) + 1;
                % save indexes of outgoing links
                if isempty(nodeparam) || isempty(nodeparam{ind}) % reinstantiate if not a cache
                    nodeparam{ind}{r} = struct();
                end
                nodeparam{ind}{r}.weights = zeros(1,self.sn.nnodes);
                nodeparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
                for c=1:size(node.output.outputStrategy{1, r}{3},2)
                    destination = node.output.outputStrategy{1, r}{3}{c}{1};
                    weight = node.output.outputStrategy{1, r}{3}{c}{2};
                    nodeparam{ind}{r}.weights(destination.index) = weight;
                end
            case RoutingStrategy.RROBIN
                nvars(ind,R+r) = nvars(ind,R+r) + 1;
                % save indexes of outgoing links
                if isempty(nodeparam) || isempty(nodeparam{ind}) % reinstantiate if not a cache
                    nodeparam{ind}{r} = struct();
                end
                nodeparam{ind}{r}.outlinks = find(self.sn.connmatrix(ind,:));
        end
    end
end

if ~isempty(self.sn) %&& isprop(self.sn,'nvars')
    self.sn.nvars = nvars;
    self.sn.nodeparam = nodeparam;
end
end
