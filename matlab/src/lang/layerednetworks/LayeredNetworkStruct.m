function lqn=LayeredNetworkStruct()
% Data structure representation for a LayeredNetwork object
%
% Copyright (c) 2012-2023, Imperial College London
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
lqn.think = {};
lqn.sched = {};
lqn.schedid = [];
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
lqn.replacement  = [];
lqn.itemproc = {};
lqn.setuptime = {};
lqn.delayofftime = {};
lqn.calltype = sparse([]);
lqn.callpair = [];
lqn.callproc = {};
lqn.callnames = {};
lqn.callhashnames = {};
%lqn.callshortnames = {};
lqn.actpretype = sparse([]);
lqn.actposttype = sparse([]);

lqn.graph = sparse([]);
lqn.taskgraph = sparse([]);
lqn.replygraph = [];

lqn.iscache = sparse(logical([]));
lqn.iscaller = sparse([]);
lqn.issynccaller = sparse([]);
lqn.isasynccaller = sparse([]);
lqn.isref = sparse(logical([]));
lqn.isfunction = sparse(logical([]));
end
