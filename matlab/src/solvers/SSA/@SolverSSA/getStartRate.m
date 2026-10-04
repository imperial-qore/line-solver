function StartN = getStartRate(self, options)
% STARTN = GETSTARTRATE()
%
% (nstations x nclasses) rate at which a class-r job BEGINS or RESUMES holding
% a server at station i, estimated over the simulated path exactly as the
% throughput is.
%
% At a lossless station with no in-service abandonment
%
%     getStartRate == getAvgTput + getPreemptRate
%
% up to simulation error; the CTMC accessor of the same name reports the exact
% value, so the two are compared with a two-sample t-test rather than an
% equality.
%
% See also getPreemptRate, SolverCTMC/getStartRate.
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
    line_error(mfilename,'This solver run produced no START rates; they are available for lang=''matlab'' runs only.');
end
StartN = self.result.startRate;
end
