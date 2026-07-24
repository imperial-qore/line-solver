function [capacity, classcap, droprule] = refreshCapacity(self)
% [CAPACITY, CLASSCAP, DROPRULE] = REFRESHCAPACITY()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
%I = getNumberOfStatefulNodes(self);
M = getNumberOfStations(self);
K = getNumberOfClasses(self);
C = self.sn.nchains;
% set zero buffers for classes that are disabled
classcap = Inf*ones(M,K);
chaincap = Inf*ones(M,K);
capacity = zeros(M,1);
droprule = DropStrategy.WAITQ*ones(M,K); 
sn = self.sn;
njobs = sn.njobs;
rates = sn.rates;
% A chain holding a spawn-target class is not population-conserving -- see _kb/04-networkstruct.md
spawnFedChain = false(1, C);
for r0 = 1:K
    if isprop(self.classes{r0}, 'spawnClass') && ~isempty(self.classes{r0}.spawnClass)
        tIdx = self.classes{r0}.spawnClass.index;
        for c0 = 1:C
            if any(sn.inchain{c0} == tIdx)
                spawnFedChain(c0) = true;
            end
        end
    end
end
for c = 1:C
    inchain = sn.inchain{c};
    for r = inchain
        chainCap = sum(njobs(inchain));
        if spawnFedChain(c)
            chainCap = Inf;
        end
        for ist=1:M
            station = self.getStationByIndex(ist);
            
            if sn.nodetype(sn.stationToNode(ist)) ~= NodeType.Source
                % Guard against dropRule array not having entry for class r; a 0 gap-fill is not a DropStrategy member (WAITQ is -1)
                if length(station.dropRule) >= r && ~isempty(station.dropRule(r)) ...
                        && station.dropRule(r) ~= 0
                    stationDropRule = station.dropRule(r);
                else
                    stationDropRule = [];
                end
                % Explicit setDropRule(WAITQ) at a finite buffer for an open class is rejected (no solver honours it) -- see _kb/04-networkstruct.md
                isUserRuleNode = ~isa(station, 'Place') && ~isa(station, 'Join');
                % classCap(r)==0 is "class absent from station", not a real finite buffer -- must not trip the WAITQ gate below
                classCapIsFinite = length(station.classCap) >= r ...
                    && station.classCap(r) > 0 && ~isinf(station.classCap(r));
                stationCapIsFinite = ~isempty(station.cap) && station.cap >= 0 ...
                    && ~isinf(station.cap) && station.cap < intmax;
                if isUserRuleNode && ~isempty(stationDropRule) && stationDropRule == DropStrategy.WAITQ ...
                        && isinf(njobs(r)) && (stationCapIsFinite || classCapIsFinite)
                    line_error(mfilename, sprintf(['Station ''%s'' declares setDropRule(WAITQ) for the open class ''%s'' at a finite capacity. ' ...
                        'LINE does not implement waiting-room blocking for an open arrival at a plain finite buffer: no solver honours this combination ' ...
                        '(SolverCTMC drops the arrival, SolverJMT blocks at the Source and ignores the capacity). ' ...
                        'Use setDropRule(DropStrategy.DROP) for a loss station (M/M/1/K), or one of the blocking policies LINE implements ' ...
                        '(DropStrategy.BAS, DropStrategy.BBS, DropStrategy.RSRD) for blocking between stations.'], ...
                        station.getName(), self.classes{r}.getName()));
                end
                % Derived from the final station capacity so it doesn't depend on setCapacity/setService call order
                if isempty(stationDropRule)
                    if isinf(station.cap) || station.cap >= intmax
                        % No buffer: the rule is never consulted.
                        stationDropRule = DropStrategy.WAITQ;
                    elseif station.cap >= 0 && ~isinf(njobs(r))
                        % Finite buffer + CLOSED class: WAITQ (the correct default is unsettled across LINE's solvers) -- see _kb/04-networkstruct.md
                        stationDropRule = DropStrategy.WAITQ;
                    else
                        % Finite buffer + OPEN class (or cap<0, JMT's unlimited sentinel): DROP, matching the CTMC reference -- see _kb/04-networkstruct.md
                        stationDropRule = DropStrategy.DROP;
                    end
                end
                droprule(ist,r) = stationDropRule;
            end
            if isnan(rates(ist,r)) && sn.nodetype(sn.stationToNode(ist)) ~= NodeType.Place
                classcap(ist,r) = 0;
                chaincap(ist,c) = 0;
            else
                chaincap(ist,c) = chainCap;
                classcap(ist,r) = chainCap;
                % Guard against classCap array not having entry for class r
                if length(station.classCap) >= r && station.classCap(r) >= 0
                    classcap(ist,r) = min(classcap(ist,r), station.classCap(r));
                end
                if station.cap >= 0
                    classcap(ist,r) = min(classcap(ist,r), station.cap);
                end
                % Finite orbit bounds population at (servers + orbit capacity): a retrial station has no waiting room -- see _kb/04-networkstruct.md
                if isprop(station,'orbitMaxJobs') && length(station.orbitMaxJobs) >= r ...
                        && station.orbitMaxJobs(r) >= 0
                    nsrv = sn.nservers(ist);
                    if ~isfinite(nsrv)
                        nsrv = 1;
                    end
                    classcap(ist,r) = min(classcap(ist,r), nsrv + station.orbitMaxJobs(r));
                end
            end
        end
    end
end
for ist=1:M
    station = self.getStationByIndex(ist);
    % If station has explicit finite cap set, use it directly (Kendall K notation)
    % Otherwise use minimum of chain cap sum and class cap sum
    if station.cap >= 0 && ~isinf(station.cap)
        % Explicit capacity set - use directly as total capacity
        capacity(ist,1) = station.cap;
    else
        capacity(ist,1) = min([sum(chaincap(ist,:)),sum(classcap(ist,:))]);
    end
end
self.sn.cap = capacity;
self.sn.classcap = classcap;
self.sn.droprule = droprule;
end
