function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = solver_ssa_nrm_space_analyzer(sn, options)
% SOLVER_SSA_ANALYZER_NRM  Performance indices from SSA/NRM simulation
%   This variant runs the next‑reaction–method SSA on the network SN and
%   returns mean throughput (XN), utilisation (UN), queue length (QN),
%   response time (RN), throughput per station (TN) and class residence
%   time (CN).  Results are based on the empirical state probabilities
%   obtained from SOLVER_SSA_NRM.
%   *** Only stations whose scheduling policy is INF, EXT or PS are
%   supported. Any other policy triggers an error. Cache nodes are not
%   supported in this simplified variant. ***

% Validate scheduling policies ---------------------------------------------------
allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS];
if any(~arrayfun(@(s) any(s == allowedSched), sn.sched))
    error('solver_ssa_analyzer_nrm:UnsupportedPolicy', ...
        'Only SchedStrategy.INF, SchedStrategy.EXT and SchedStrategy.PS are supported.');
end

line_debug('SSA NRM space analyzer starting: nstations=%d, nclasses=%d', sn.nstations, sn.nclasses);

end
