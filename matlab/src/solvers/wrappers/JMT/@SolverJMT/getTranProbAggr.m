function [Pi_t, SSnode_a] = getTranProbAggr(self, node)
% [PI_T, SSNODE_A] = GETTRANPROBAGGR(NODE)
%
% Empirical transient probability of per-class aggregate states at a station,
% estimated from iter_max independent JMT replications.
%
% Output:
%   Pi_t     - [t, pi_t] where pi_t(k,j) = empirical fraction of replications
%              whose aggregate state at the requested node at time t(k) equals
%              SSnode_a(j,:).
%   SSnode_a - Unique aggregate state vectors observed at NODE across all
%              replications and time points (sorted, one row per state, one
%              column per class).
if GlobalConstants.DummyMode
    Pi_t = NaN;
    SSnode_a = NaN;
    return
end
options = self.getOptions;
initSeed = self.options.seed;
if isfield(options,'timespan')  && isfinite(options.timespan(2))
    sn = self.getStruct;
    ist = sn.nodeToStation(node.index);
    if ist < 1
        line_error(mfilename,'getTranProbAggr in SolverJMT requires the input node to be a station.');
    end
    if sn.nodetype(node.index) == NodeType.Source
        line_error(mfilename,'getTranProbAggr in SolverJMT does not apply to Source nodes.');
    end

    tu = [];
    TranSysStateAggr = cell(options.iter_max,1);
    for it=1:options.iter_max
        self.options.seed = initSeed + it - 1;
        TranSysStateAggr{it} = self.sampleSysAggr;
        tu = union(tu, TranSysStateAggr{it}.t);
    end
    self.options.seed = initSeed; % restore in case of interruption
    tu = tu(:);
    nT = length(tu);
    iterMax = options.iter_max;
    nclasses = sn.nclasses;

    % interpolate per-station state onto union grid for the requested station
    nodeStates = NaN(nT, nclasses, iterMax);
    for it=1:iterMax
        tit = TranSysStateAggr{it}.t;
        block = TranSysStateAggr{it}.state{ist};
        if isempty(block) || isempty(tit) || any(isinf(block(:)))
            continue
        end
        nodeStates(:,:,it) = interp1(tit, block, tu, 'previous');
    end

    % flatten to (nT*iterMax) x nclasses, drop rows with any NaN, take unique
    flat = reshape(permute(nodeStates,[1,3,2]), nT*iterMax, nclasses);
    validMask = ~any(isnan(flat),2);
    SSnode_a = unique(flat(validMask,:), 'rows');
    nStates = size(SSnode_a,1);

    % empirical fraction at each time point
    pi_t = zeros(nT, nStates);
    if nStates > 0
        for it=1:iterMax
            block = nodeStates(:,:,it);
            [~, idx] = ismember(block, SSnode_a, 'rows');
            for t_idx = 1:nT
                if idx(t_idx) > 0
                    pi_t(t_idx, idx(t_idx)) = pi_t(t_idx, idx(t_idx)) + 1;
                end
            end
        end
        pi_t = pi_t / iterMax;
    end

    Pi_t = [tu, pi_t];
else
    line_error(mfilename,'getTranProbAggr in SolverJMT requires to specify a finite timespan T, e.g., SolverJMT(model,''timespan'',[0,T]).');
end
end
