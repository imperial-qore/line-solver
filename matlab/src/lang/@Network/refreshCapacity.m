function [capacity, classcap, dropid] = refreshCapacity(self)
% [CAPACITY, CLASSCAP, DROPRULE] = REFRESHCAPACITY()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
%I = getNumberOfStatefulNodes(self);
M = getNumberOfStations(self);
K = getNumberOfClasses(self);
C = self.sn.nchains;
% set zero buffers for classes that are disabled
classcap = Inf*ones(M,K);
chaincap = Inf*ones(M,K);
capacity = zeros(M,1);
dropid = -ones(M,K); %DropStrategy.WaitingBuffer
sn = self.sn;
njobs = sn.njobs;
rates = sn.rates;
for c = 1:C
    inchain = sn.inchain{c};
    for r = inchain
        chainCap = sum(njobs(inchain));
        for i=1:M
            station = self.getStationByIndex(i);
            
            if sn.nodetype(sn.stationToNode(i)) ~= NodeType.ID_SOURCE     
                dropid(i,r) = station.dropRule(r);
            end 
            if isnan(rates(i,r)) && sn.nodetype(sn.stationToNode(i)) ~= NodeType.ID_PLACE
                classcap(i,r) = 0;
                chaincap(i,c) = 0;
            else
                chaincap(i,c) = chainCap;
                classcap(i,r) = chainCap;
                if station.classCap(r) >= 0
                    classcap(i,r) = min(classcap(i,r), station.classCap(r));
                end
                if station.cap >= 0
                    classcap(i,r) = min(classcap(i,r), station.cap);
                end
            end
        end
    end
end
for i=1:M
    capacity(i,1) = min([sum(chaincap(i,:)),sum(classcap(i,:))]);
end
self.sn.cap = capacity;
self.sn.classcap = classcap;
self.sn.dropid = dropid;
end
