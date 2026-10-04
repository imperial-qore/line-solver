function session = ag_cluster_session(exec, agents, reset)
% SESSION = AG_CLUSTER_SESSION(EXEC, AGENTS, RESET)
%
% Open (or reuse) the connections to the ag-worker processes and assign the
% agents to them once, up front.
%
% WHY THE ASSIGNMENT IS SEPARATE FROM THE SWEEP. An agent's static data -- its
% local rate matrix and the passive/active matrices of the actions it takes
% part in -- does not change across the fixed point; only the reversed rates x
% do, and those are one double per action. Shipping the matrices once and then
% exchanging x per sweep is what makes remote execution worth doing; shipping
% them every sweep would put an O(N^2) payload on the wire to save an O(N^3)
% solve and lose most of the benefit.
%
% The partition is round-robin over the worker list in agent index order, so it
% is a pure function of (number of agents, number of workers): a rerun assigns
% the same agents to the same workers, which is what keeps a distributed run
% reproducible.
%
% A worker that cannot be reached is NOT fatal. Its slot is left empty and
% ag_cluster_sweep solves those agents locally, because any agent can be solved
% anywhere given x. A cluster run therefore degrades to a slower run, never to
% a wrong one or a failed one.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

persistent CACHE

if nargin >= 3 && reset
    if ~isempty(CACHE)
        closeAll(CACHE);
    end
    CACHE = [];
    session = [];
    return;
end

key = sprintf('%s|', exec.endpoints{:});
key = sprintf('%s#%d', key, numel(agents));
if ~isempty(CACHE) && strcmp(CACHE.key, key)
    session = CACHE;
    return;
end
if ~isempty(CACHE)
    closeAll(CACHE);
end

nw = numel(exec.endpoints);
session = struct('key', key, 'sock', {cell(1, nw)}, 'in', {cell(1, nw)}, ...
    'out', {cell(1, nw)}, 'owns', {cell(1, nw)}, 'live', false(1, nw), ...
    'timeout', exec.timeout);

% Round-robin partition, computed before any connection is attempted so that a
% dead worker does not shift the others' agents.
for k = 1:numel(agents)
    w = mod(k - 1, nw) + 1;
    session.owns{w}(end+1) = k;
end

for w = 1:nw
    ep = exec.endpoints{w};
    parts = strsplit(ep, ':');
    if numel(parts) ~= 2
        line_warning(mfilename, 'Malformed AG worker endpoint ''%s'', expected host:port. Its agents run locally.\n', ep);
        continue;
    end
    host = parts{1};
    port = str2double(parts{2});
    try
        sock = java.net.Socket();
        sock.connect(java.net.InetSocketAddress(host, port), round(session.timeout * 1000));
        sock.setSoTimeout(round(session.timeout * 1000));
        out = java.io.PrintWriter(java.io.OutputStreamWriter(sock.getOutputStream(), 'UTF-8'), true);
        in = java.io.BufferedReader(java.io.InputStreamReader(sock.getInputStream(), 'UTF-8'));

        msg = struct('op', 'assign', 'agents', {agentPayload(agents, session.owns{w})});
        out.println(jsonencode(msg));
        reply = readReply(in);
        if ~isfield(reply, 'op') || ~strcmp(reply.op, 'assigned')
            error('ag:assign', 'worker did not acknowledge the assignment');
        end

        session.sock{w} = sock;
        session.in{w} = in;
        session.out{w} = out;
        session.live(w) = true;
    catch ME
        line_warning(mfilename, ['AG worker %s is unreachable (%s). Its %d agent(s) run ' ...
            'locally instead.\n'], ep, ME.message, numel(session.owns{w}));
    end
end

CACHE = session;

end

function payload = agentPayload(agents, ks)
% The static half of each owned agent, in the wire form the worker expects.
payload = cell(1, numel(ks));
for t = 1:numel(ks)
    payload{t} = agents(ks(t));
end
end

function reply = readReply(in)
line = in.readLine();
if isempty(line)
    error('ag:eof', 'worker closed the connection');
end
reply = jsondecode(char(line));
end

function closeAll(session)
for w = 1:numel(session.sock)
    if ~isempty(session.sock{w})
        try
            session.out{w}.println(jsonencode(struct('op', 'bye')));
            session.sock{w}.close();
        catch
            % A worker that already went away needs no goodbye.
        end
    end
end
end
