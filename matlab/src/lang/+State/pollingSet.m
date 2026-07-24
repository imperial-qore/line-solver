function space_var = pollingSet(pinfo, space_var, pos, swk, ctr)
% SPACE_VAR = POLLINGSET(PINFO, SPACE_VAR, POS, SWK, CTR)
%
% Write the polling controller into the local-variable columns of a state row.
% Columns that State.pollingInfo elides are dropped: they are reconstructible
% from the rest of the state (see State.pollingGet), so materializing them
% would split each state into copies that no observation can tell apart.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if pinfo.ipos > 0
    space_var(:, pinfo.ipos) = pos;
end
if pinfo.iswk > 0
    space_var(:, pinfo.iswk) = swk;
end
if pinfo.ictr > 0
    space_var(:, pinfo.ictr) = ctr;
end
end
