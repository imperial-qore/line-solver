function Prob = ldesHistProb(space, time, R, stations, targets)
% PROB = LDESHISTPROB(SPACE, TIME, R, STATIONS, TARGETS)
% Residence-time probability of an aggregate joint state in an LDES histogram.
%
% SPACE is the (nstates x nstations*R) joint-state matrix and TIME the (nstates
% x 1) residence time of each row, as ldesHistogram returns them; R is
% sn.nclasses, which fixes the column stride. STATIONS is a vector of 1-based
% station indices to constrain and TARGETS a cell array of matching per-class
% job-count vectors; a target shorter than R constrains only the classes it
% names. Returns 0 when the state is never visited, which is a measurement.
%
% The layout is the station-major, class-minor one ctmc_state_space_aggr builds:
% column (i-1)*R+k is the number of class-k jobs at station i.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Prob = 0;
if isempty(space) || isempty(time)
    return
end
total = sum(time(:));
if total <= 0
    return
end

match = true(size(space, 1), 1);
for s = 1:numel(stations)
    ist = stations(s);
    tgt = targets{s}(:).';
    L = min(numel(tgt), R);
    cols = (ist - 1) * R + (1:L);
    if any(cols > size(space, 2))
        line_error(mfilename, sprintf('SolverLDES: station %d is past the end of the state histogram, which holds %d columns for %d classes.', ist, size(space, 2), R));
    end
    match = match & all(abs(space(:, cols) - repmat(tgt(1:L), size(space, 1), 1)) < 1e-9, 2);
end

Prob = sum(time(match)) / total;
end
