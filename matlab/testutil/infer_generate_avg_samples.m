function [respT, arvR, util] = infer_generate_avg_samples(model, samples, C)
%generateAvgSamples Use SSA solver to generate average response, arrival
%and utilization data for given network
    numStations = size(model.getNodes(), 1)
    respT = zeros(samples, numStations, C);
    arvR = zeros(samples, numStations, C);
    util = zeros(samples, numStations);

    for s=1:samples
        solver = SolverSSA(model, 'samples', 100000, 'seed', s);
        respT(s, :, :) = solver.getAvgRespT();
        arvR(s, :, :) = solver.getAvgArvR();
        util(s, :) = sum(solver.getAvgUtil(), 2);
    end    
end

