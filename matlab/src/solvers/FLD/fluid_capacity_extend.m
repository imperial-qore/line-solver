function con = fluid_capacity_extend(con, stg, terms)
% CON = FLUID_CAPACITY_EXTEND(CON, STG, TERMS)
%
% Extend every cap to the staging coordinates that hold mass INSIDE it.
%
% A waiting room is outside the region it feeds, which is the whole point of it
% -- but it is not outside every OTHER limit. Where two regions overlap, an
% admission into the inner one is an INTERNAL move of the outer one: the job
% leaves a station of the outer region, waits, and re-enters a station of the
% same outer region, never having left it. Counting only the state coordinates
% would take that mass out of the outer cap for as long as it waits, so the outer
% cap would be met on paper while the region actually held more; and the Newton
% system that results is inconsistent rather than merely inexact -- two
% overlapping regions stalled at residual 5e-1 with the inner cap exceeded.
%
% A room counts toward a row when the row weighs the room's DESTINATION and also
% weighs every station that FEEDS it -- that is exactly "the job was inside and
% stays inside". A room fed from outside is a queue at the door and counts
% nowhere, as before.
%
% Parameters:
%   con   - constraint set from FLUID_CAPACITY_CONSTRAINTS
%   stg   - staging from FLUID_CAPACITY_STAGING
%   terms - the event representation from FLUID_MOMENT_TERMS
%
% Returns:
%   con - the same struct with As (ncon x stg.n) filled in
%
% See also FLUID_CAPACITY_CONSTRAINTS, FLUID_CAPACITY_STAGING.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ncon = numel(con.b);
con.As = zeros(ncon, stg.n);
if ncon == 0 || stg.n == 0
    return
end
tol = sqrt(GlobalConstants.Zero);
Dn = min(terms.D, 0);
Dp = max(terms.D, 0);

% which stations feed each room, and which coordinate its mass lands on
feeds = cell(stg.n,1);
dest = zeros(stg.n,1);
for e = find(stg.adm(:))'
    j = stg.admStage(e);
    if j < 1
        continue
    end
    feeds{j} = union(feeds{j}, find(Dn(:,e) < -tol)');
    landing = find(Dp(:,e) > tol)';
    for s = landing
        if con.member(stg.region(j), s)
            dest(j) = s;
            break
        end
    end
end

for c = 1:ncon
    for j = 1:stg.n
        if dest(j) < 1 || isempty(feeds{j})
            continue
        end
        w = con.A(c, dest(j));
        if w <= 0
            continue
        end
        if all(con.A(c, feeds{j}) > 0)
            con.As(c,j) = w;
        end
    end
end
end
