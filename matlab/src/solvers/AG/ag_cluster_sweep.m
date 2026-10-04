function [pi, Q] = ag_cluster_sweep(x, Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta)
% [PI, Q] = AG_CLUSTER_SWEEP(X, AA, PB, L, ACT, PSV, NUMPROCESSES, A, N, META)
%
% One sweep of the reversed-rate fixed point with the agents solved on remote
% ag-worker processes.
%
% THE ONLY THING THAT CROSSES THE WIRE PER SWEEP IS X, one double per action,
% and back come the agents' stationary vectors. The generator itself is rebuilt
% here rather than shipped, because assembling it is the cheap half (O(N^2))
% and solving it is the expensive half (O(N^3)); see ag_agent_generator.
%
% A worker that is missing, slow or broken is not fatal: its agents are solved
% locally through the same agent path, so the result is the run's result either
% way and only the wall clock changes. That is a property of the decomposition,
% not a fallback bolted on -- an agent depends on the rest of the model only
% through x, so anyone holding x can solve it.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

persistent AGENTS

Q = cell(1, numProcesses);
pi = cell(1, numProcesses);

if isempty(AGENTS) || numel(AGENTS) ~= numProcesses
    AGENTS = ag_cluster_agents(Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta);
end
session = ag_cluster_session(meta.exec, AGENTS);

% Every agent's generator is rebuilt here in any case: the metrics stage reads
% Q, and rebuilding is cheaper than transporting it.
for k = 1:numProcesses
    Q{k} = ag_agent_generator(k, x, Aa, Pb, L, ACT, PSV, A, N);
end

pending = true(1, numProcesses);

for w = 1:numel(session.live)
    ks = session.owns{w};
    if isempty(ks) || ~session.live(w)
        continue;
    end
    try
        msg = struct('op', 'sweep', 'mode', 'finite', 'x', x(:)');
        session.out{w}.println(jsonencode(msg));
        line = session.in{w}.readLine();
        if isempty(line)
            error('ag:eof', 'worker closed the connection mid-sweep');
        end
        reply = jsondecode(char(line));
        if ~isfield(reply, 'op') || ~strcmp(reply.op, 'swept')
            error('ag:proto', 'unexpected reply ''%s''', reply.op);
        end
        for t = 1:numel(reply.agents)
            r = reply.agents(t);
            pi{r.k} = r.pi(:)';
            pending(r.k) = false;
        end
    catch ME
        line_warning(mfilename, ['AG worker %d failed mid-sweep (%s); its %d agent(s) ' ...
            'are solved locally for the rest of the run.\n'], w, ME.message, numel(ks));
        session.live(w) = false;
    end
end

% Whatever no worker answered for, solve here. This covers an unreachable
% worker, a worker lost mid-run and the agents of a malformed endpoint alike.
for k = find(pending)
    pi{k} = ag_solve_component(Q{k}, meta.mph(k), meta.nlev(k), meta.level{k});
end

end
