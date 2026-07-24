function [arvR, util, avgQLengths, timeIntervals] = infer_generate_transient_samples(model, samples, C, stride)
%generateTransientAvgSamples Use SSA solver to generate steady-state arrival
%rates and transient utilization and queue-length data for given network

    solver = SolverSSA(model);

    [~,~,~,~] = solver.getAvg();
    arvR = repmat(solver.getAvgArvR(), samples, 1);

    [Q, U, T] = model.getTranHandles();

    [QNt,Ut,~] = SolverSSA(model,'force', true, 'timespan',[0,5]).getTranAvg(Q,U,T);

    stationCount = model.getNumberOfNodes();
    jobCount = size(model.classes, 1);

    avgQLengths = cell(stationCount, jobCount);
    util = cell(stationCount, jobCount);

    timeIntervals = QNt{1,1}.t;
    timeIntervals = timeIntervals(1:stride:size(timeIntervals,1));

    for i=1:stationCount
        for c=1:jobCount
            met = QNt{i,c}.metric;
            len = size(met, 1);
            avgQLengths{i, c} = met(1:stride:len);
            umet = Ut{i,c}.metric;
            util{i, c} = umet(1:stride:size(umet,1));
        end
    end
end

