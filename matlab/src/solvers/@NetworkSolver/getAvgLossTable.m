function LossTable = getAvgLossTable(self)
% LOSSTABLE = GETAVGLOSSTABLE()
%
% Table of loss (drop) metrics for every station-class pair that receives
% offered traffic. For each such pair it reports the offered arrival rate
% (ArvR), the carried throughput (Tput), the loss rate (ArvR - Tput, i.e. the
% rate of jobs dropped by finite capacity, blocking, or reneging) and the loss
% ratio (LossRate / ArvR, the fraction of offered jobs that are lost).
%
% Only station-class pairs with ArvR > 0 are listed: the loss identity
% LossRate = ArvR - Tput holds where a class is actually offered to a station,
% and it excludes the Source, whose offered arrival rate is zero. A lossless
% station has ArvR = Tput and therefore LossRate = LossRatio = 0.
%
% See also NetworkSolver.getAvgTable, NetworkSolver.getAvg
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[~,~,~,TN,AN] = self.getAvg();
sn = self.model.getStruct();

% Fork-Join quorum sibling-drop rate (LDES only), station-indexed. At a
% synchronizing Join the identity LossRate = ArvR - Tput does not hold (Tput is
% in parent units, discarded siblings in sibling units), so on Join rows the
% explicit drop rate replaces ArvR - Tput. ArvR is already the offered sibling
% rate (flow balance over all forked siblings), so LossRatio = drop / ArvR.
DropRateJoin = [];
if isprop(self, 'result') && isstruct(self.result) && isfield(self.result, 'DropRateJoin')
    DropRateJoin = self.result.DropRateJoin;
end

Station = {};
JobClass = {};
ArvR = [];
Tput = [];
LossRate = [];
LossRatio = [];
for ist = 1:size(AN,1)
    for r = 1:size(AN,2)
        a = AN(ist,r);
        if ~(isfinite(a) && a > 0)
            continue
        end
        t = TN(ist,r);
        d = 0;
        if ~isempty(DropRateJoin) && ist <= size(DropRateJoin,1) && r <= size(DropRateJoin,2)
            d = DropRateJoin(ist,r);
        end
        if isfinite(d) && d > 0
            lr = d;
            lc = d / a;
        else
            lr = a - t;
            lc = (a - t) / a;
        end
        Station{end+1,1} = sn.nodenames{sn.stationToNode(ist)}; %#ok<AGROW>
        JobClass{end+1,1} = sn.classnames{r}; %#ok<AGROW>
        ArvR(end+1,1) = a; %#ok<AGROW>
        Tput(end+1,1) = t; %#ok<AGROW>
        LossRate(end+1,1) = lr; %#ok<AGROW>
        LossRatio(end+1,1) = lc; %#ok<AGROW>
    end
end

LossTable = Table(Station, JobClass, ArvR, Tput, LossRate, LossRatio);
LossTable = IndexedTable(LossTable);
end
