function ldes_require_finite_pas_buffers(model)
% LDES_REQUIRE_FINITE_PAS_BUFFERS(MODEL)
% Refuse an OI/PAS station with no finite buffer before it reaches LDES.
%
% mu(c) cannot cross JSON as a function, so linemodel_save materializes it as a
% macrostate table, and LDES needs a finite buffer to bound the state whether or
% not the rate saturates: saturation only keeps the table compact, it cannot
% bound the queue length. Without this gate an uncapped station was written over
% the 10-per-open-class box and served clamped at the box edge without a word.
% The saturation search is the shared pas_saturation_cutoffs (matlab/src/io).
%
% The reference refusal is Python save_model's (linemodel_io.py); C++ raises the
% same two messages in io::require_finite_pas_buffers, called by
% ldes::solver_ldes. The JAR needs no twin: its NetworkSolver constructor builds
% the struct WITH the initial state, so FromMarginal already refuses an uncapped
% PAS station for every solver. The messages are the reference's word for word,
% naming the station. The gate is at SOLVE time, NOT in linemodel_save, which
% must still save such a model. model.getStruct() with the initial state refuses
% the same model generically (State.fromMarginalAndStarted, 'PAS stations require
% finite capacity ...', not naming the station), so the gate has to run before
% it to be heard: NetworkSolver.getAvg (getAvg is Sealed, hence a SolverLDES
% branch there), SolverLDES.getStruct (getAvgTable), SolverLDES.runAnalyzer
% (getTranAvg) and SolverLDES.solveCli (every engine run, including the entry
% points that build the struct without the initial state). getProb*, sample*,
% getAvgReward, getTranReward and initFromSolver take their struct from
% SolverLDES.getStruct, never from model.getStruct, so that one gate covers
% them; a new SolverLDES method must do the same. A buffer is finite when
% sn.cap is, so a closed station bounded by its chain population passes.
% see _kb/09-ldes-and-cache.md
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isa(model, 'Network')
    return
end
sn = model.getStruct(false);
K = sn.nclasses;
nodes = model.getNodes;
for ind = 1:numel(nodes)
    node = nodes{ind};
    if ~isa(node, 'Queue') || isempty(node.svcRateFun)
        continue
    end
    ist = sn.nodeToStation(ind);
    if isnan(ist) || ist < 1 || isfinite(sn.cap(ist))
        continue
    end
    [tau, bounded] = pas_saturation_cutoffs(node.svcRateFun, K);
    if bounded
        tauStr = ['[', strjoin(arrayfun(@(x) sprintf('%d', x), tau, 'UniformOutput', false), ', '), ']'];
        line_error(mfilename, sprintf(['PAS station ''%s'' has no finite buffer: LDES state-space ' ...
            'generation requires one. Its order-independent service rate saturates at class ' ...
            'counts %s, so call setCap(N) with N comfortably above the mean occupancy.'], ...
            node.name, tauStr));
    end
    line_error(mfilename, sprintf(['PAS station ''%s'' has an unbounded order-independent service ' ...
        'rate (mu does not saturate) and no finite buffer; call setCap(N) on the station.'], node.name));
end
end

