function [pos, swk, ctr] = pollingGet(pinfo, space_var, srvclass)
% [POS, SWK, CTR] = POLLINGGET(PINFO, SPACE_VAR, SRVCLASS)
%
% Read the polling controller out of the local-variable columns of a single
% state row. SRVCLASS is the class currently in the service facility, or 0 when
% the facility is empty.
%
% The columns that State.pollingInfo elides are reconstructed here, so that
% callers always see a complete controller:
%   pos  when not materialized, every switchover is immediate, so the server is
%        either serving (and then it stands at the buffer of the job in
%        service) or parked (and then its position is unobservable: from a park
%        with only immediate switchovers the walk reaches any buffer in zero
%        time, so all positions have identical dynamics and 1 is canonical).
%   swk  when not materialized, no switchover takes time, so the server is
%        never inside one.
%   ctr  when not materialized the discipline is EXHAUSTIVE, which bounds a
%        visit by the buffer draining rather than by a budget.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if pinfo.ipos > 0
    pos = space_var(1, pinfo.ipos);
elseif srvclass > 0
    pos = srvclass;
else
    pos = 1;
end
if pinfo.iswk > 0
    swk = space_var(1, pinfo.iswk);
else
    swk = 0;
end
if pinfo.ictr > 0
    ctr = space_var(1, pinfo.ictr);
else
    ctr = 0;
end
end
