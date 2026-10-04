function StartN = getStartRate(self, options)
% STARTN = GETSTARTRATE()
%
% (nstations x nclasses) rate at which a class-r job BEGINS or RESUMES holding
% a server at station i, i.e. pi*F*e over the START filtration.
%
% At a lossless station with no in-service abandonment
%
%     getStartRate == getAvgTput + getPreemptRate
%
% because every job starts service once per entry into a server and every
% preemption is followed by exactly one later resume or restart. At a
% non-preemptive station this collapses to startRate == throughput.
%
% This is an accessor on SolverCTMC, not a MetricType: it adds no column to
% getAvgTable.
%
% See also getPreemptRate, getEventFiltration.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = self.getOptions;
end
if isempty(self.result) || ~isfield(self.result,'startRate') || isempty(self.result.startRate)
    self.runAnalyzer(options);
end
if ~isfield(self.result,'startRate') || isempty(self.result.startRate)
    line_error(mfilename,'This solver run produced no START rates; they are unavailable for lang=''cpp'' results taken from line-cli.');
end
StartN = self.result.startRate;
end
