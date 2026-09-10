function [rows, probs] = pollingLand(pinfo, q, mode, budget, space_buf, space_srv, space_var, K, Ks, pieist, R)
% [ROWS, PROBS] = POLLINGLAND(PINFO, Q, MODE, BUDGET, SPACE_BUF, SPACE_SRV, SPACE_VAR, K, KS, PIEIST, R)
%
% Materialize the state rows a polling server lands in after State.pollingNext
% has resolved (Q, MODE, BUDGET), together with the probability of each. The
% three inputs SPACE_BUF/SPACE_SRV/SPACE_VAR are single rows describing the
% station at the instant the decision is taken, i.e. with the completed job (if
% any) already removed from the service facility.
%
% PROBS splits a landing across the entry phases of a phase-type: which phase a
% service or a switchover starts in is a random choice, so one decision yields
% one row per entry phase, weighted by the corresponding entry probability.
% Callers fold PROBS into the transition rate rather than into outprob, since
% the branching happens at the instant the active event fires.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

rows = [];
probs = [];
switch mode
    case 1 % start or continue a visit at q: pull a waiting class-q job into service
        buf = space_buf;
        buf(1,q) = buf(1,q) - 1;
        pentry = pieist{q};
        for kentry = 1:K(q)
            if pentry(kentry) <= 0
                continue
            end
            srv = space_srv;
            srv(1,Ks(q)+kentry) = srv(1,Ks(q)+kentry) + 1;
            var = State.pollingSet(pinfo, space_var, q, 0, budget);
            rows(end+1,:) = [buf, srv, var]; %#ok<AGROW>
            probs(end+1,1) = pentry(kentry); %#ok<AGROW>
        end
    case 2 % enter the switchover into q: the facility stays empty while walking
        swpie = pinfo.swpie{q};
        for kentry = 1:pinfo.Ksw(q)
            if swpie(kentry) <= 0
                continue
            end
            var = State.pollingSet(pinfo, space_var, q, kentry, 0);
            rows(end+1,:) = [space_buf, space_srv, var]; %#ok<AGROW>
            probs(end+1,1) = swpie(kentry); %#ok<AGROW>
        end
    case 0 % park: held until the next arrival, see State.pollingNext
        var = State.pollingSet(pinfo, space_var, q, 0, 0);
        rows = [space_buf, space_srv, var];
        probs = 1;
end
end
