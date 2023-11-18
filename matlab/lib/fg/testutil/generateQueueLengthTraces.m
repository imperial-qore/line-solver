function [timeIntervals, queueLengthMatrix] = generateQueueLengthTraces(model, stride)
%GENERATEQUEUELENGTHTRACES Use JMT solver to generate queue length traces
%for given network

    solver = SolverJMT(model);

    [~,~,~,~] = solver.getAvg();

    [Qt,Ut,Tt] = model.getTranHandles();

    [QNt,~,~] = SolverJMT(model, 'force', true, 'timespan',[0,5]).getTranAvg(Qt,Ut,Tt);
    
    stationCount = model.getNumberOfNodes();
    jobCount = size(model.classes, 1);

    avgQLengths = cell(stationCount, jobCount);

    timeIntervals = QNt{1,1}.t;
    timeIntervals = timeIntervals(1:stride:size(timeIntervals,1));

    for i=1:stationCount
        for c=1:jobCount
            met = QNt{i,c}.metric;
            len = size(met, 1);
            avgQLengths{i, c} = met(1:stride:len);
        end
    end

    queueLengthMatrix = avgQLengths;
end

