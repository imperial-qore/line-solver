classdef ModelAdapter < Copyable
    % Static class to transform and adapt models, providing functionality for:
    % - Creating tagged job models for response time analysis
    % - Fork-join network transformations (formerly from api/fj/)
    % - Model preprocessing and adaptation operations
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)

        model = PMIF2LINE(filename,modelName);
        sn = PMIF2QN(filename,verbose);

        % ========== Model adaptation ==========
        [taggedModel, taggedJob] = tagChain(model, chain, jobclass, suffix);
        newmodel = removeClass(model, jobclass);
        [chainModel, alpha, deaggInfo] = aggregateChains(model, suffix);
        [fesModel, fesStation, deaggInfo] = aggregateFES(model, stationSubset, options);

        % ========== Fork-Join Methods (formerly from api/fj/) ==========
        ri = findPaths(sn, P, start, endNode, r, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel);
        ri = findPathsCS(sn, P, curNode, endNode, curClass, toMerge, QN, TN, currentTime, fjclassmap, fjforkmap, nonfjmodel, visited);        
        ri = paths(sn, P, start, endNode, toMerge, QN, TN, currentTime);
        [ri, stat, RN] = pathsCS(sn, orignodes, P, f, joinIdx, r, RN, currentTime, toMerge);
        [forks, parents] = sortForks(sn, nonfjstruct, fjforkmap, fjclassmap, nonfjmodel);
        [nonfjmodel, fjclassmap, fjforkmap, fanout, prov] = mmt(model, forkLambda);
        ok = refreshServicesFromBase(nonfjmodel, prov);
        [nonfjmodel, fjclassmap, fjforkmap, fj_auxiliary_delays] = ht(model);
        [fjmodel, fjsn, fjclassmap, fjforkmap, fjjoinmap, fjbranchmap, fjtagmap] = fjtag(model);
    end

    % The converters live on the path in matlab/src/io, not in this class folder, so a bare
    % declaration left ModelAdapter.X undefined; each static method forwards to the path function.
    methods (Static)
        function model = JMT2LINE(varargin)
            model = JMT2LINE(varargin{:});
        end
        function model = JMVA2LINE(varargin)
            model = JMVA2LINE(varargin{:});
        end
        function model = JSIM2LINE(varargin)
            model = JSIM2LINE(varargin{:});
        end
        function java_model = LINE2JLINE(line_model)
            java_model = LINE2JLINE(line_model);
        end
        function model = QN2LINE(varargin)
            model = QN2LINE(varargin{:});
        end
        function lqnmodel = QN2LQN(varargin)
            lqnmodel = QN2LQN(varargin{:});
        end

        % ========== Model I/O ==========
        function LINE2MATLAB(varargin)
            LINE2MATLAB(varargin{:});
        end
        function LINE2JAVA(varargin)
            LINE2JAVA(varargin{:});
        end
        function QN2MATLAB(varargin)
            QN2MATLAB(varargin{:});
        end
        function QN2JAVA(varargin)
            QN2JAVA(varargin{:});
        end
        function varargout = LQN2MATLAB(varargin)
            % generator form LQN2MATLAB(lqnmodel, modelName, fid), or the legacy loader form
            [varargout{1:nargout}] = LQN2MATLAB(varargin{:});
        end
        function LQN2JAVA(varargin)
            LQN2JAVA(varargin{:});
        end
    end

end
