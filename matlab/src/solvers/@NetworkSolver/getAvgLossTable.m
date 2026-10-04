function varargout = getAvgLossTable(self,varargin)
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
% A Join is the exception, and its row is computed differently: see
% sn_join_droprate. A standard join reports LossRate = 0 there; a quorum join
% reports the rate at which stragglers are discarded.
%
% See also NetworkSolver.getAvgTable, NetworkSolver.getAvg
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% The result recorder captures the returned table together with the solver
% that produced it, so cross-codebase parity is asserted against the values a
% solver RETURNED rather than the text it printed. Off unless a run asked for
% it (LineResultRecorder.enable), and then it costs one appdata lookup here.
% The wrapper exists so that recording happens on EVERY exit path, including
% the early returns inside the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getAvgLossTable_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'loss', varargout{1});
end

function LossTable = getAvgLossTable_impl(self)
% GETAVGLOSSTABLE_IMPL Implementation of GETAVGLOSSTABLE; see the wrapper above.

[~,~,~,TN,AN] = self.getAvg();
sn = self.model.getStruct();

% Fork-Join sibling-drop rate, station-indexed. At a synchronizing Join the
% identity LossRate = ArvR - Tput does not hold, because the two rates are in
% different units: ArvR counts the SIBLINGS offered (N per parent job) and Tput
% the PARENT jobs released. Reading ArvR - Tput there charges (N-1)/N of the
% offered traffic as lost at EVERY join, standard joins included. On a Join row
% the drop rate therefore replaces ArvR - Tput unconditionally: the solver's own
% measurement when it supplies one (SolverLDES counts the discards on its sample
% path), otherwise sn_join_droprate's ArvR - K*Tput, which is exact given the
% two rates. LossRatio stays drop / ArvR.
DropRateJoin = [];
if isprop(self, 'result') && isstruct(self.result) && isfield(self.result, 'DropRateJoin')
    DropRateJoin = self.result.DropRateJoin;
end
isJoinRow = false(size(AN,1),1);
for ist0 = 1:size(AN,1)
    ind0 = sn.stationToNode(ist0);
    isJoinRow(ist0) = ind0 >= 1 && sn.nodetype(ind0) == NodeType.Join;
end
if any(isJoinRow) && (isempty(DropRateJoin) || ~isequal(size(DropRateJoin), size(AN)))
    DropRateJoin = sn_join_droprate(sn, TN, AN);
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
        if isJoinRow(ist)
            lr = max(0, d);
            lc = lr / a;
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
