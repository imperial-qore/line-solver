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
for c = 1:C
    inchain = sn.inchain{c};
    for r = inchain
        chainCap = sum(njobs(inchain));
        for ist=1:M
            station = self.getStationByIndex(ist);
            
            if sn.nodetype(sn.stationToNode(ist)) ~= NodeType.Source
                % Guard against dropRule array not having entry for class r
                if length(station.dropRule) >= r && ~isempty(station.dropRule(r))
                    stationDropRule = station.dropRule(r);
                else
                    stationDropRule = [];
                end
                % Default to Drop for finite capacity queues if no explicit drop rule set
                if isempty(stationDropRule)
                    if ~isinf(station.cap) && station.cap < intmax
                        stationDropRule = DropStrategy.DROP;
                    else
                        stationDropRule = DropStrategy.WAITQ;
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
