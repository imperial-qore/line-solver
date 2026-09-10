function agents = ag_cluster_agents(Aa, Pb, L, ACT, PSV, numProcesses, A, N, meta)
% AGENTS = AG_CLUSTER_AGENTS(AA, PB, L, ACT, PSV, NUMPROCESSES, A, N, META)
%
% The STATIC description of every agent, i.e. everything a worker needs that
% does not change across the fixed point. Only the reversed rates x change, so
% this is built once per solve and shipped once per worker.
%
% Matrices travel as [row, col, value] triplets with 1-based indices, because
% the passive and active matrices of an action are nearly empty (an action
% touches one level transition) and a dense N-by-N payload would be mostly
% zeros. The receiving worker rebuilds the same matrices and therefore the same
% generator, which is what makes a remote agent's answer identical to a local
% one rather than merely close.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

agents = struct('k', {}, 'N', {}, 'mph', {}, 'nlev', {}, 'level', {}, ...
    'L', {}, 'passive', {}, 'active', {});

for k = 1:numProcesses
    a = struct();
    a.k = k;
    a.N = N(k);
    a.mph = meta.mph(k);
    a.nlev = meta.nlev(k);
    a.level = meta.level{k}(:)';
    a.L = triplets(L{k});

    pas = struct('c', {}, 'M', {});
    act = struct('c', {}, 'M', {});
    for c = 1:A
        if PSV(c) == k
            pas(end+1) = struct('c', c, 'M', {triplets(Pb{c})}); %#ok<AGROW>
        elseif ACT(c) == k
            act(end+1) = struct('c', c, 'M', {triplets(Aa{c})}); %#ok<AGROW>
        end
    end
    a.passive = pas;
    a.active = act;

    agents(end+1) = a; %#ok<AGROW>
end

end

function t = triplets(M)
% Non-zero entries of M as an n-by-3 [row col value] array, 1-based. An empty
% matrix yields a 0-by-3, which jsonencode writes as [] and the worker reads
% back as "no entries" rather than as a malformed payload.
if isempty(M)
    t = zeros(0, 3);
    return;
end
[i, j, v] = find(M);
t = [i(:), j(:), v(:)];
end
