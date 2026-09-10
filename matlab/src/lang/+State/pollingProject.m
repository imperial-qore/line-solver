function block = pollingProject(pinfo, trips)
% BLOCK = POLLINGPROJECT(PINFO, TRIPS)
%
% Project full [pos, swk, ctr] controller triples onto the columns that
% State.pollingInfo materializes. The elided columns are reconstructible from
% the rest of the state (see State.pollingGet), so keeping them would split
% each state into copies no observation can tell apart.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

keep = logical([pinfo.wpos, pinfo.wswk, pinfo.wctr]);
block = trips(:, keep);
end
