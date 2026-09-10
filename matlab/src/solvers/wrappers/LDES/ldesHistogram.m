function [space, time] = ldesHistogram(solver)
% [SPACE, TIME] = LDESHISTOGRAM(SOLVER)
% Exact joint-state residence-time histogram of one LDES run.
%
% Runs the fully JSON-mediated solve with --export-histogram and returns the
% (nstates x nstations*nclasses) joint-state matrix and the (nstates x 1)
% residence time of each row, so P(state) = t(state)/sum(t) is exact on the
% sampled path.
%
% WHY THIS AND NOT THE TRAJECTORY. The transient QNt series the engine returns
% under --trajectory holds INTERVAL TIME-AVERAGES of the queue length, not the
% integer states the path visits, so it cannot be compared against a state to
% estimate a probability. Every probability the LDES wrapper reports goes
% through this histogram, as getAvgReward already does.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

space = [];
time = [];
data = solver.solveCli(solver.getOptions, {'--export-histogram'});
if isstruct(data) && isfield(data, 'stateHistogram') && ~isempty(data.stateHistogram)
    space = ldesJson2mat(ldesGetField(data.stateHistogram, 'space', []), [], []);
    time = ldesJson2mat(ldesGetField(data.stateHistogram, 'time', []), [], []);
end
if isempty(space) || isempty(time)
    line_error(mfilename, 'SolverLDES: the engine returned no state histogram, so no state probability can be read from this run.');
end
end
