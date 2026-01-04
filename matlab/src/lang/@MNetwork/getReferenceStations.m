function refstat = getReferenceStations(self)
% REFSTAT = GETREFERENCESTATIONS()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

K = getNumberOfClasses(self);
refstat = zeros(K,1);
classes = self.classes;
for k=1:K    
    refstat(k,1) = classes{k}.refstat.stationIndex;
    %refstat(k,1) =  findstring(getStationNames(self),self.classes{k}.refstat.name);
end
end
