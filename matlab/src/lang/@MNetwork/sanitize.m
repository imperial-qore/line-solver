function sanitize(self)
% SANITIZE()

% Preprocess model to ensure consistent parameterization.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Convert Signal placeholders to OpenSignal or ClosedSignal based on routing
self.convertSignalPlaceholders();

% Initialize Cache node properties always (even if checks disabled)
% This is required for proper cache functionality
K = getNumberOfClasses(self);
for ist=1:self.getNumberOfNodes
    if isa(self.nodes{ist}, 'Cache')
        for k=1:K
            if k > length(self.nodes{ist}.popularity) || isempty(self.nodes{ist}.popularity{k})
                self.nodes{ist}.popularity{k} = Disabled.getInstance();
            end
        end
        if isempty(self.nodes{ist}.accessProb)
            self.nodes{ist}.accessProb = cell(K,self.nodes{ist}.items.nitems);
            for v=1:K
                for k=1:self.nodes{ist}.items.nitems
                    % accessProb{v,k}(l,p) is the cost (probability) for a user-v request to item k in list l to access list p
                    if isempty(self.nodes{ist}.graph)
                        self.nodes{ist}.accessProb{v,k} = diag(ones(1,self.nodes{ist}.nLevels),1);
                        self.nodes{ist}.accessProb{v,k}(1+self.nodes{ist}.nLevels,1+self.nodes{ist}.nLevels) = 1;
                    else
                        self.nodes{ist}.accessProb{v,k} = self.nodes{ist}.graph{k};
                    end
                end
            end
        end
        self.nodes{ist}.server.hitClass = round(self.nodes{ist}.server.hitClass);
        self.nodes{ist}.server.missClass = round(self.nodes{ist}.server.missClass);
    end
end

% Skip remaining validation if checks are disabled
if ~self.enableChecks
    return;
end

if isempty(self.sn)
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    for ist=1:self.getNumberOfNodes
        switch class(self.nodes{ist})
            case 'Logger'
                %no-op
            case 'ClassSwitch'
                %no-op
            case 'Join'
                for k=1:K
                    self.nodes{ist}.classCap(k) = Inf;
                    self.nodes{ist}.dropRule(k) = DropStrategy.WAITQ;
                end
            case 'Queue'
                hasEnabledService = false;
                for k=1:K
                    if k > length(self.nodes{ist}.server.serviceProcess) || isempty(self.nodes{ist}.server.serviceProcess{k})
                        self.nodes{ist}.serviceProcess{k} = Disabled.getInstance();
                        self.nodes{ist}.classCap(k) = 0;
                        self.nodes{ist}.dropRule(k) = DropStrategy.WAITQ;
                        self.nodes{ist}.schedStrategyPar(k) = 0;
                        self.nodes{ist}.server.serviceProcess{k} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                    else
                        % Check if this class has an enabled service
                        if ~isa(self.nodes{ist}.server.serviceProcess{k}{3}, 'Disabled')
                            hasEnabledService = true;
                        end
                    end
                end
                % Also check for heterogeneous services
                if ~hasEnabledService && isprop(self.nodes{ist}, 'heteroServiceDistributions') && ~isempty(self.nodes{ist}.heteroServiceDistributions)
                    if self.nodes{ist}.heteroServiceDistributions.Count > 0
                        hasEnabledService = true;
                    end
                end
                % Validate that at least one class has service configured
                % Skip validation for models with Cache, Petri net, or Fork/Join elements
                hasSpecialElements = false;
                for jnode=1:self.getNumberOfNodes
                    if isa(self.nodes{jnode}, 'Cache') || isa(self.nodes{jnode}, 'Place') || isa(self.nodes{jnode}, 'Transition') || isa(self.nodes{jnode}, 'Fork') || isa(self.nodes{jnode}, 'Join')
                        hasSpecialElements = true;
                        break;
                    end
                end
                if ~hasEnabledService && ~hasSpecialElements
                    line_error(mfilename, sprintf('Queue ''%s'' has no service configured for any job class. Use setService() to configure service times.', self.nodes{ist}.name));
                end
                for k=1:K
                    if isa(self.nodes{ist}.serviceProcess{k},'Disabled')
                        self.nodes{ist}.setRouting(self.classes{k},RoutingStrategy.DISABLED);
                    end
                end                
                switch SchedStrategy.toId(self.nodes{ist}.schedStrategy)
                    case SchedStrategy.SEPT
                        servTime = zeros(1,K);
                        for k=1:K
                            servTime(k) = self.nodes{ist}.serviceProcess{k}.getMean;
                        end
                        if length(unique(servTime)) ~= K
                            line_error(mfilename,'SEPT does not support identical service time means.');
                        end
                        [servTimeSorted] = sort(unique(servTime));
                        self.nodes{ist}.schedStrategyPar = zeros(1,K);
                        for k=1:K
                            if ~isnan(servTime(k))
                                self.nodes{ist}.schedStrategyPar(k) = find(servTimeSorted == servTime(k));
                            else
                                self.nodes{ist}.schedStrategyPar(k) = find(isnan(servTimeSorted));
                            end
                        end
                    case SchedStrategy.LEPT
                        servTime = zeros(1,K);
                        for k=1:K
                            servTime(k) = self.nodes{ist}.serviceProcess{k}.getMean;
                        end
                        if length(unique(servTime)) ~= K
                            line_error(mfilename,'LEPT does not support identical service time means.');
                        end
                        [servTimeSorted] = sort(unique(servTime),'descend');
                        self.nodes{ist}.schedStrategyPar = zeros(1,K);
                        for k=1:K
                            self.nodes{ist}.schedStrategyPar(k) = find(servTimeSorted == servTime(k));
                        end
                end
            case 'Delay'
                hasEnabledService = false;
                for k=1:K
                    if k > length(self.nodes{ist}.server.serviceProcess) || isempty(self.nodes{ist}.server.serviceProcess{k})
                        self.nodes{ist}.serviceProcess{k} = Disabled.getInstance();
                        self.nodes{ist}.classCap(k) = 0;
                        self.nodes{ist}.dropRule(k) = DropStrategy.WAITQ;
                        self.nodes{ist}.schedStrategyPar(k) = 0;
                        self.nodes{ist}.server.serviceProcess{k} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                    else
                        % Check if this class has an enabled service
                        if ~isa(self.nodes{ist}.server.serviceProcess{k}{3}, 'Disabled')
                            hasEnabledService = true;
                        end
                    end
                end
                % Validate that at least one class has service configured
                % Skip validation for models with Cache, Petri net, or Fork/Join elements
                hasSpecialElements = false;
                for jnode=1:self.getNumberOfNodes
                    if isa(self.nodes{jnode}, 'Cache') || isa(self.nodes{jnode}, 'Place') || isa(self.nodes{jnode}, 'Transition') || isa(self.nodes{jnode}, 'Fork') || isa(self.nodes{jnode}, 'Join')
                        hasSpecialElements = true;
                        break;
                    end
                end
                if ~hasEnabledService && ~hasSpecialElements
                    line_error(mfilename, sprintf('Delay ''%s'' has no service configured for any job class. Use setService() to configure service times.', self.nodes{ist}.name));
                end
                for k=1:K
                    if isa(self.nodes{ist}.serviceProcess{k},'Disabled')
                        self.nodes{ist}.setRouting(self.classes{k},RoutingStrategy.DISABLED);
                    end
                end                
                switch SchedStrategy.toId(self.nodes{ist}.schedStrategy)
                    case SchedStrategy.SEPT
                        servTime = zeros(1,K);
                        for k=1:K
                            servTime(k) = self.nodes{ist}.serviceProcess{k}.getMean;
                        end
                        [~,self.nodes{ist}.schedStrategyPar] = sort(servTime);
                end
            case 'Sink'
                for k=1:K
                    self.getSink.setRouting(self.classes{k},RoutingStrategy.DISABLED);
                end
            case 'Source'
                for k=1:K
                    if k > length(self.nodes{ist}.input.sourceClasses) || isempty(self.nodes{ist}.input.sourceClasses{k})
                        self.nodes{ist}.input.sourceClasses{k} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                        self.nodes{ist}.setArrival(self.classes{k},Disabled.getInstance());                        
                    end
                end
            case 'Place'
                for k=1:K
                    if k > length(self.nodes{ist}.server.serviceProcess) || isempty(self.nodes{ist}.server.serviceProcess{k})
                        self.nodes{ist}.schedStrategyPar(k) = 0;
                        self.nodes{ist}.server.serviceProcess{k} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                    end
                end
            case 'Transition'
                for k=1:K
                    if k > length(self.nodes{ist}.server.serviceProcess) || isempty(self.nodes{ist}.server.serviceProcess{k})
                        %                         self.nodes{i}.schedStrategyPar(k) = 0;
                        self.nodes{ist}.server.serviceProcess{k} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                    end
                end
        end
    end
    for ist=1:M
        if isempty(self.getIndexSourceStation) || ist ~= self.getIndexSourceStation
            for r=1:K
                switch self.stations{ist}.server.className
                    case 'ServiceTunnel'
                        % do nothing
                    case 'Cache'
                        self.stations{ist}.setProbRouting(self.classes{r}, self.stations{ist}, 0.0);
                    otherwise
                        if isempty(self.stations{ist}.server.serviceProcess{r})
                            self.stations{ist}.server.serviceProcess{r} = {[],ServiceStrategy.LI,Disabled.getInstance()};
                        end
                end
            end
        end
    end

    % Check if model has Cache, Petri net, or Fork/Join elements
    hasSpecialElements = false;
    for inode=1:self.getNumberOfNodes
        if isa(self.nodes{inode}, 'Cache') || isa(self.nodes{inode}, 'Place') || isa(self.nodes{inode}, 'Transition') || isa(self.nodes{inode}, 'Fork') || isa(self.nodes{inode}, 'Join')
            hasSpecialElements = true;
            break;
        end
    end

    % Validate that each job class has service configured at least one station
    % Skip validation for models with Cache, Petri net, or Fork/Join elements
    if ~hasSpecialElements
        for k=1:K
            hasServiceAtAnyStation = false;
            isUsedInCacheSetRead = false;
            isOpenClass = isa(self.classes{k}, 'OpenClass');

            % Check if job class is used in setRead at a Cache node
            for inode=1:self.getNumberOfNodes
                if isa(self.nodes{inode}, 'Cache')
                    if size(self.nodes{inode}.popularity,2) >= k && ~isempty(self.nodes{inode}.popularity{k}) && ~isa(self.nodes{inode}.popularity{k}, 'Disabled')
                        isUsedInCacheSetRead = true;
                        break;
                    end
                end
            end

            for ist=1:M
                if k <= length(self.stations{ist}.server.serviceProcess) && ...
                   ~isempty(self.stations{ist}.server.serviceProcess{k}) && ...
                   length(self.stations{ist}.server.serviceProcess{k}) >= 3 && ...
                   ~isa(self.stations{ist}.server.serviceProcess{k}{3}, 'Disabled')
                    hasServiceAtAnyStation = true;
                    break;
                end
            end

            % Only closed classes need explicit setService configuration (open classes route to Sink by design)
            if ~hasServiceAtAnyStation && ~isUsedInCacheSetRead && ~isOpenClass
                line_error(mfilename, sprintf('Job class ''%s'' has no service configured at any station. Every job class must have service configured at least one station using setService().', self.classes{k}.name));
            end
        end
    end
end
end
