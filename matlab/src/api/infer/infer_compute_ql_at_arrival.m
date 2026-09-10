function ql = infer_compute_ql_at_arrival(at, at_jobid, rt, rt_jobid, class, R)
% INFER_COMPUTE_QL_AT_ARRIVAL Compute per-class queue lengths at arrival.
%
% Reconstructs the queue state seen by each arriving job using arrival
% and departure times. Arrival times and response times are matched by
% job ID, so they need not be in the same order or come from the same
% data source.
%
% At ties, departures are processed before arrivals.
%
% Inputs:
%   at        - arrival times (column vector, n x 1)
%   at_jobid  - job IDs for arrival times (column vector, n x 1)
%   rt        - response times (column vector, m x 1, m >= n)
%   rt_jobid  - job IDs for response times (column vector, m x 1)
%   class     - class of each arrival sample (column vector, n x 1,
%               ordered consistently with at and at_jobid)
%   R         - number of classes
%
% Returns:
%   ql        - n x R matrix of per-class queue lengths at each arrival,
%               rows ordered consistently with the input at/at_jobid
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

n = length(at);

% Match response times to arrivals by job ID
[found, loc] = ismember(at_jobid, rt_jobid);
if ~all(found)
    error('infer_compute_ql_at_arrival: not all arrival job IDs found in response time job IDs.');
end
rt_matched = rt(loc);

% Sort arrivals by time
[at_sorted, sortIdx] = sort(at);
class_sorted = class(sortIdx);
rt_sorted = rt_matched(sortIdx);

exitTimes = at_sorted + rt_sorted;

% Event list: [time, type (-1=dep/+1=arv), sortedIdx, class]
events = zeros(2*n, 4);
events(1:n, :) = [at_sorted, ones(n,1), (1:n)', class_sorted];
events(n+1:2*n, :) = [exitTimes, -ones(n,1), (1:n)', class_sorted];

% Sort by time; departures (-1) before arrivals (+1) at same time
events = sortrows(events, [1, 2]);

state = zeros(1, R);
ql_sorted = zeros(n, R);
for i = 1:2*n
    c = events(i, 4);
    if events(i, 2) == 1  % arrival
        state(c) = state(c) + 1;
        ql_sorted(events(i, 3), :) = state;
    else  % departure
        state(c) = state(c) - 1;
    end
end

% Unsort back to original input order
ql = zeros(n, R);
ql(sortIdx, :) = ql_sorted;

end
