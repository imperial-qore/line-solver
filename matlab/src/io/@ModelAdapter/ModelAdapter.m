classdef ModelAdapter < Copyable
    % Static class to transform and adapt models, providing functionality for:
    % - Creating tagged job models for response time analysis
    % - Fork-join network transformations (formerly from api/fj/)
    % - Model preprocessing and adaptation operations
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)

        model = JMT2LINE(filename,modelName);
        model = JMVA2LINE(filename,modelName);
        model = JSIM2LINE(filename,modelName);
        java_model = LINE2JLINE(line_model);
        model = PMIF2LINE(filename,modelName);
        sn = PMIF2QN(filename,verbose);
        model = QN2LINE(sn, modelName);
        lqnmodel=QN2LQN(model);

        % ========== Model I/O ==========
        LINE2MATLAB(model, filename);
        LINE2JAVA(model, filename);
        QN2MATLAB(model, modelName, fid);
        QN2JAVA(model, modelName, fid, headers);
        LQN2MATLAB(lqnmodel, modelName, fid);
        LQN2JAVA(lqnmodel, modelName, fid);

        % ========== Model adaptation ==========
        [taggedModel, taggedJob] = tagChain(model, chain, jobclass, suffix);
        newmodel = removeClass(model, jobclass);
        [chainModel, alpha, deaggInfo] = aggregateChains(model, suffix);
        [fesModel, fesStation, deaggInfo] = aggregateFES(model, stationSubset, options);

        % ========== Fork-Join Methods (formerly from api/fj/) ==========
        ri = findPaths(sn, P, start, endNode, r, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel);
        ri = findPathsCS(sn, P, curNode, endNode, curClass, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel);        
        ri = paths(sn, P, start, endNode, toMerge, QN, TN, currentTime);
        [ri, stat, RN] = pathsCS(sn, orignodes, P, f, joinIdx, r, RN, currentTime, toMerge);
        [forks, parents] = sortForks(sn, nonfjstruct, fjforkmap, fjclassmap, nonfjmodel);
        [nonfjmodel, fjclassmap, fjforkmap, fanout] = mmt(model, forkLambda);
        [nonfjmodel, fjclassmap, fjforkmap, fj_auxiliary_delays] = ht(model);
    end

end