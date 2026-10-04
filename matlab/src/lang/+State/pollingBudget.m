function budget = pollingBudget(pinfo, nbufq)
% BUDGET = POLLINGBUDGET(PINFO, NBUFQ)
%
% Initial value of the ctr column of a visit that starts at a buffer holding
% NBUFQ waiting jobs. NBUFQ is the class-q population at the polling instant:
% the service facility is empty when a visit starts, so it is the whole
% class-q population at the station.
%
% See State.pollingInfo for the meaning of ctr under each discipline.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch pinfo.ptype
    case PollingType.EXHAUSTIVE
        budget = 0; % unused: the visit ends when the buffer drains
    case PollingType.GATED
        budget = nbufq; % serve exactly the jobs found at the polling instant
    case PollingType.KLIMITED
        budget = pinfo.pk; % serve at most K, fewer if the buffer drains first
    case PollingType.DECREMENTING
        budget = nbufq - 1; % serve until the population drops one below the level found
    otherwise
        line_error(mfilename, sprintf('Unsupported polling type: %d.', pinfo.ptype));
end
end
