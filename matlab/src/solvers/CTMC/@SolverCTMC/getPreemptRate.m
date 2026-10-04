function PreemptN = getPreemptRate(self, options)
% PREEMPTN = GETPREEMPTRATE()
%
% (nstations x nclasses) rate at which a class-r job HOLDING A SERVER at
% station i is pushed back into the buffer, i.e. pi*F*e over the PREEMPT
% filtration. It is identically zero at a non-preemptive station.
%
% Preempt-resume and preempt-independent stations report the SAME rate: which
% phase the displaced job resumes in is not a property of how often it is
% displaced.
%
% This is an accessor on SolverCTMC, not a MetricType: it adds no column to
% getAvgTable.
%
% See also getStartRate, getEventFiltration.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = self.getOptions;
end
if isempty(self.result) || ~isfield(self.result,'preemptRate') || isempty(self.result.preemptRate)
    self.runAnalyzer(options);
end
if ~isfield(self.result,'preemptRate') || isempty(self.result.preemptRate)
    line_error(mfilename,'This solver run produced no PREEMPT rates; they are unavailable for lang=''cpp'' results taken from line-cli.');
end
PreemptN = self.result.preemptRate;
end
