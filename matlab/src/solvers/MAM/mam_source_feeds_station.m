function tf = mam_source_feeds_station(sn, jst, ist)
% TF = MAM_SOURCE_FEEDS_STATION(SN, JST, IST)
%
% True when station JST is a Source whose whole output IS the whole arrival
% stream of station IST.
%
% The exact GI/M/c closed forms (QSYS_PHMC, QSYS_DMC) take an INTER-ARRIVAL law
% and read it out of some station's sn.proc. That substitution is only legitimate
% when the arrivals at IST are, sample path for sample path, the epochs that
% station generates: JST must be a Source (so its process is an inter-arrival
% law and not a service law), all of its flow must reach IST directly, and
% nothing else may reach IST. Anywhere else in a network the input of a station
% is a DEPARTURE process, which is not the upstream station's service law and
% which only the decomposition produces.
%
% Both tests are written over SN_RT_STATIONS rather than over sn.rt, because
% sn.rt is indexed by stateful node and a Router, Cache or stateful class switch
% between the Source and IST would otherwise read as a foreign feeder.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tf = false;
M = sn.nstations;
K = sn.nclasses;
if jst < 1 || jst > M || ist < 1 || ist > M || jst == ist
    return;
end
if sn.nodetype(sn.stationToNode(jst)) ~= NodeType.Source
    return;
end

rtst = sn_rt_stations(sn);
rows = @(st) ((st-1)*K + (1:K));
srcRows = rows(jst);
dstCols = rows(ist);

% All of the source's flow reaches IST directly.
intoDst = sum(rtst(srcRows, dstCols), 2);
total = sum(rtst(srcRows, :), 2);
if sum(total) <= GlobalConstants.FineTol
    return;
end
if any(abs(intoDst - total) > GlobalConstants.FineTol)
    return;
end

% Nothing else reaches IST.
for kst = 1:M
    if kst == jst
        continue;
    end
    if sum(sum(rtst(rows(kst), dstCols))) > GlobalConstants.FineTol
        return;
    end
end

tf = true;
end
