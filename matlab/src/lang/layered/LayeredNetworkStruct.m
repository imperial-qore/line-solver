function lqn=LayeredNetworkStruct()
% Data structure representation for a LayeredNetwork object
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lqn=[]; %faster than lqn=struct();
lqn.nidx = 0;  % total number of hosts, tasks, entries, and activities, except the reference tasks
lqn.nhosts = 0;
lqn.ntasks = 0;
lqn.nentries = 0;
lqn.nacts = 0;
lqn.ncalls = 0;
lqn.hshift = 0;
lqn.tshift = 0;
lqn.eshift = 0;
lqn.ashift = 0;
lqn.cshift = 0;
lqn.nidx = 0;
lqn.tasksof = {};
lqn.entriesof = {};
lqn.actsof = {};
lqn.callsof = {};
lqn.hostdem = {};
lqn.hostdem_type = [];      % ProcessType enum vector
lqn.hostdem_params = {};    % Cell array of parameter vectors
lqn.hostdem_mean = [];      % Precomputed mean vector
lqn.hostdem_scv = [];       % Precomputed SCV vector
lqn.hostdem_proc = {};      % MAP/PH process {D0,D1} matrices

lqn.actthink = {};
lqn.actthink_type = [];
lqn.actthink_params = {};
lqn.actthink_mean = [];
lqn.actthink_scv = [];
lqn.actthink_proc = {};

lqn.think = {};
lqn.think_type = [];
lqn.think_params = {};
lqn.think_mean = [];
lqn.think_scv = [];
lqn.think_proc = {};
lqn.sched = [];
lqn.names = {};
lqn.hashnames = {};
%lqn.shortnames = {};
lqn.mult = [];
lqn.repl = [];
lqn.type = [];
%lqn.replies = [];
lqn.parent = [];

lqn.nitems = [];
lqn.itemcap  = {};
lqn.replacestrat  = [];

lqn.itemproc = {};
lqn.itemproc_type = [];
lqn.itemproc_params = {};
lqn.itemproc_mean = [];
lqn.itemproc_scv = [];
lqn.itemproc_proc = {};

lqn.setuptime = {};
lqn.setuptime_type = [];
lqn.setuptime_params = {};
lqn.setuptime_mean = [];
lqn.setuptime_scv = [];
lqn.setuptime_proc = {};

lqn.delayofftime = {};
lqn.delayofftime_type = [];
lqn.delayofftime_params = {};
lqn.delayofftime_mean = [];
lqn.delayofftime_scv = [];
lqn.delayofftime_proc = {};
lqn.calltype = sparse([]);
lqn.callpair = [];

lqn.callproc = {};
lqn.callproc_type = [];
lqn.callproc_params = {};
lqn.callproc_mean = [];
lqn.callproc_scv = [];
lqn.callproc_proc = {};

lqn.callnames = {};
lqn.callhashnames = {};

% Open arrival distributions (populated in getStruct for entries with arrivals)
lqn.arrival = {};
lqn.arrival_type = [];
lqn.arrival_params = {};
lqn.arrival_mean = [];
lqn.arrival_scv = [];
lqn.arrival_proc = {};
%lqn.callshortnames = {};
lqn.actpretype = sparse([]);
lqn.actposttype = sparse([]);

lqn.graph = sparse([]);
lqn.taskgraph = sparse([]);
lqn.replygraph = [];
lqn.actphase = [];  % Phase number (1 or 2) for each activity

lqn.iscache = sparse(logical([]));
lqn.iscaller = sparse([]);
lqn.issynccaller = sparse([]);
lqn.isasynccaller = sparse([]);
lqn.isref = sparse(logical([]));
lqn.isfunction = sparse(logical([]));
end
