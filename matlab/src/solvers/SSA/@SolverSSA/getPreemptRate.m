function PreemptN = getPreemptRate(self, options)
% PREEMPTN = GETPREEMPTRATE()
%
% (nstations x nclasses) rate at which a class-r job HOLDING A SERVER at
% station i is pushed back into the buffer, estimated over the simulated path.
% It is identically zero at a non-preemptive station.
%
% See also getStartRate, SolverCTMC/getPreemptRate.
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
    line_error(mfilename,'This solver run produced no PREEMPT rates; they are available for lang=''matlab'' runs only.');
end
PreemptN = self.result.preemptRate;
end
