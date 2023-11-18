function [respT, arvR, util] = generateTransientAvgSamples(model, samples, C)
%generateAvgSamples Use JMT solver to generate average response, arrival
%and utilization data for given network
    
    solver = SolverJMT(model);

    [~,~,~,~] = solver.getAvg();
    arvR = repmat(solver.getAvgArvR(), samples, 1);


    [Q, U, T] = model.getTranHandles();

    [QNt,Ut,Tt] = SolverJMT(model,'force', true, 'timespan',[0,5]).getTranAvg(Q,U,T);
    
    stationCount = model.getNumberOfNodes();
    jobCount = size(model.classes, 1);

    avgQLengths = cell(stationCount, jobCount);

    timeIntervals = QNt{1,1}.t;
    timeIntervals = timeIntervals(1:stride:size(timeIntervals,1));

    for i=1:stationCount
        util = Ut{i,c}.metric;
        for c=1:jobCount
            met = QNt{i,c}.metric;
            len = size(met, 1);
            avgQLengths{i, c} = met(1:stride:len);
        end

    end
end

