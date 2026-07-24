function [timeIntervals, queueLengthMatrix] = infer_generate_qlen_traces(model, stride)
%GENERATEQUEUELENGTHTRACES Use SSA solver to generate queue length traces
%for given network

    sn = model.getStruct();
    nStateful = sn.nstateful;
    jobCount = sn.nclasses;

    solver = SolverSSA(model, 'force', true, 'timespan', [0, Inf], ...
        'samples', 10000, 'method', 'serial');
    ts = solver.sampleSysAggr(10000);

    % Subsample by stride
    nSamples = length(ts.t);
    idx = 1:stride:nSamples;
    timeIntervals = ts.t(idx);

    avgQLengths = cell(nStateful, jobCount);
    for isf = 1:nStateful
        for c = 1:jobCount
            avgQLengths{isf, c} = ts.state{isf}(idx, c);
        end
    end

    queueLengthMatrix = avgQLengths;
end

