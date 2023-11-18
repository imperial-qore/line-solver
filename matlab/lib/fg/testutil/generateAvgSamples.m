function [respT, arvR, util] = generateAvgSamples(model, samples, C)
%generateAvgSamples Use JMT solver to generate average response, arrival
%and utilization data for given network
    numStations = size(model.getNodes(), 1)
    respT = zeros(samples, numStations, C);
    arvR = zeros(samples, numStations, C);
    util = zeros(samples, numStations);

    for s=1:samples
        solver = SolverJMT(model, 'samples', 100000, 'seed', s);
        respT(s, :, :) = solver.getAvgRespT();
        arvR(s, :, :) = solver.getAvgArvR();
        util(s, :) = sum(solver.getAvgUtil(), 2);
    end    
end

