% Regression test for synchronized fork-join arrivals in SolverMAM.
%
% The synchronized join arrival process is built with mmap_max. The result
% must be normalized before it is reused, otherwise mmap_lambda can report
% infeasible throughputs at the join.

sn = struct();
sn.nnodes = 5;
sn.nclasses = 1;
sn.nodetype = [NodeType.Source, NodeType.Queue, NodeType.Queue, NodeType.Join, NodeType.Sink];
sn.rtnodes = zeros(sn.nnodes * sn.nclasses);

sourceIdx = 1;
queue1Idx = 2;
queue2Idx = 3;
joinIdx = 4;
sinkIdx = 5;

sn.rtnodes(sourceIdx, queue1Idx) = 1;
sn.rtnodes(sourceIdx, queue2Idx) = 1;
sn.rtnodes(queue1Idx, joinIdx) = 1;
sn.rtnodes(queue2Idx, joinIdx) = 1;
sn.rtnodes(joinIdx, sinkIdx) = 1;

fjSyncMap = struct();
fjSyncMap.nodeSync = zeros(sn.nnodes, sn.nnodes);
fjSyncMap.nodeSync(joinIdx, queue1Idx) = 1;
fjSyncMap.nodeSync(joinIdx, queue2Idx) = 1;

config = struct();
config.fj_sync_q_len = 2;
config.merge = 'super';
config.compress = 'none';
config.space_max = 100;

DEP = cell(sn.nnodes, sn.nclasses);
DEP{sourceIdx,1} = map_exponential(20);
DEP{queue1Idx,1} = map_exponential(20);
DEP{queue2Idx,1} = map_exponential(1 / 0.049878);
DEP{joinIdx,1} = map_exponential(20);

ARV = solver_mam_traffic_mmap(sn, DEP, config, fjSyncMap);
joinArrival = ARV{joinIdx};
joinRate = mmap_lambda(joinArrival);
rowSums = sum(joinArrival{1} + joinArrival{2}, 2);

assert(all(abs(rowSums) < 1e-12), ...
    'solver_mam_traffic_mmap returned a non-feasible MMAP at the join.');
assert(joinRate(1) > 0, ...
    'solver_mam_traffic_mmap returned zero synchronized throughput at the join.');
assert(joinRate(1) < 0.05, ...
    'solver_mam_traffic_mmap returned a join throughput above the branch throughput.');
