function [filt, ev] = getEventFiltration(self, eventType, options)
% [FILT, EV] = GETEVENTFILTRATION(EVENTTYPE)
%
% Filtration of a DERIVED event type: an {nstations x nclasses} cell of sparse
% (state x state) matrices whose (s,ns) entry is the rate at which the
% transition s -> ns carries one such event at that station for that class.
% EV holds the matching Event descriptors, in the same cell layout.
%
% EVENTTYPE must be EventType.START or EventType.PREEMPT. The two are not
% synchronizations: they are tags on the ARV and DEP arcs that cause them, so
% they are NOT part of getGenerator's eventFilt (which is paired one-to-one
% with sn.sync and summed as D1 by sample.m) and are kept here instead.
%
% Example:
%   solver = SolverCTMC(model);
%   F = solver.getEventFiltration(EventType.PREEMPT);
%   pi = solver.getProbSys();       % preemption rate of class r at station i
%   rate_ir = pi * F{i,r} * ones(size(F{i,r},2),1);
%
% See also getStartRate, getPreemptRate.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3
    options = self.getOptions;
end
if nargin < 2 || isempty(eventType)
    line_error(mfilename,'getEventFiltration requires an event type, either EventType.START or EventType.PREEMPT.');
end
if eventType ~= EventType.START && eventType ~= EventType.PREEMPT
    line_error(mfilename, sprintf(['getEventFiltration serves the derived events only (START, PREEMPT); ' ...
        '%s is a synchronization and its filtration is the one getGenerator returns.'], EventType.toText(eventType)));
end

if isempty(self.result) || ~isfield(self.result,'auxFilt') || isempty(self.result.auxFilt)
    self.getGenerator(options);
end
if ~isfield(self.result,'auxFilt') || isempty(self.result.auxFilt)
    line_error(mfilename,'This model produced no derived event filtration; it is not available for lang=''cpp'' generators taken from line-cli.');
end

sn = self.getStruct;
if eventType == EventType.START
    filt = self.result.auxFilt.start;
else
    filt = self.result.auxFilt.preempt;
end

ev = cell(size(filt));
for i = 1:size(filt,1)
    ind = sn.stationToNode(i);
    for r = 1:size(filt,2)
        ev{i,r} = Event(eventType, ind, r);
    end
end
end
