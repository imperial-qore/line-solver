function varargout = getAvgRegionLossTable(self,varargin)
% LOSSTABLE = GETAVGREGIONLOSSTABLE()
%
% Table of loss (drop) metrics per finite-capacity region and class, for
% regions that drop jobs (DROP rule). Each row reports the offered arrival rate
% (carried Tput plus drop rate), the carried throughput, the loss rate (region
% drop rate) and the loss ratio (LossRate / ArvR).
%
% Only regions with offered traffic are listed. The table is empty for solvers
% that do not track region drops: only the LDES simulation populates the FCR
% drop rate (self.result.FCR.DropRateNfcr).
%
% See also NetworkSolver.getAvgLossTable, Queue.setDropRule
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
[varargout{1:max(nargout,1)}] = getAvgRegionLossTable_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'regionloss', varargout{1});
end

function LossTable = getAvgRegionLossTable_impl(self)
% GETAVGREGIONLOSSTABLE_IMPL Implementation of GETAVGREGIONLOSSTABLE; see the wrapper above.

self.getAvgTable();

Region = {};
JobClass = {};
ArvR = [];
Tput = [];
LossRate = [];
LossRatio = [];
if isprop(self, 'result') && isfield(self.result, 'FCR') ...
        && isfield(self.result.FCR, 'DropRateNfcr')
    sn = self.model.getStruct();
    TNfcr = self.result.FCR.TNfcr;
    DR = self.result.FCR.DropRateNfcr;
    for f = 1:size(DR, 1)
        for r = 1:size(DR, 2)
            t = TNfcr(f, r);
            d = DR(f, r);
            a = t + d;
            if ~(isfinite(a) && a > 0)
                continue
            end
            Region{end+1,1} = sprintf('FCRegion%d', f); %#ok<AGROW>
            JobClass{end+1,1} = sn.classnames{r}; %#ok<AGROW>
            ArvR(end+1,1) = a; %#ok<AGROW>
            Tput(end+1,1) = t; %#ok<AGROW>
            LossRate(end+1,1) = d; %#ok<AGROW>
            LossRatio(end+1,1) = d / a; %#ok<AGROW>
        end
    end
end

LossTable = Table(Region, JobClass, ArvR, Tput, LossRate, LossRatio);
LossTable = IndexedTable(LossTable);
end
