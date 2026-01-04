classdef SolverLQNS < Solver
    % A solver that interfaces the LQNS to LINE.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverLQNS(model, varargin)
            % SELF = SOLVERLQNS(MODEL, VARARGIN)
            self@Solver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if ~SolverLQNS.isAvailable() && ~self.options.config.remote
                line_error(mfilename,['SolverLQNS requires the lqns and lqsim commands to be available on the system path.\n' ...
                    'You can install them from: http://www.sce.carleton.ca/rads/lqns/\n\n' ...
                    'Alternatively, use remote execution via Docker:\n' ...
                    '  1. Pull and run: docker run -d -p 8080:8080 imperialqore/lqns-rest:latest\n' ...
                    '  2. Configure remote execution in MATLAB:\n' ...
                    '     options = SolverOptions(@()SolverLQNS);\n' ...
                    '     options.config.remote = true;\n' ...
                    '     options.config.remote_url = ''http://localhost:8080'';\n' ...
                    '     solver = SolverLQNS(model, options);']);
            end
        end

        function sn = getStruct(self)
            %GETSTRUCT Retrieve the model structure
            sn = self.model.getStruct();
        end

        function varargout = getAvg(varargin)
            %GETAVG Proxy to getEnsembleAvg
            [varargout{1:nargout}] = getEnsembleAvg(varargin{:});
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTable(self)
            % [AVGTABLE,QT,UT,RT,WT,TT] = GETAVGTABLE()
            if (GlobalConstants.DummyMode)
                [AvgTable, QT, UT, RT, TT, WT] = deal([]);
                return
            end

            if ~isempty(self.obj)
                avgTable = self.obj.getEnsembleAvg();
                [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(avgTable);
            else
                [QN,UN,RN,TN,AN,WN] = getAvg(self);
            end

            % attempt to sanitize small numerical perturbations
            variables = {QN, UN, RN, TN, AN, WN};  % Put all variables in a cell array
            for i = 1:length(variables)
                rVar = round(variables{i} * 10);
                toRound = abs(variables{i} * 10 - rVar) < GlobalConstants.CoarseTol * variables{i} * 10;
                variables{i}(toRound) = rVar(toRound) / 10;
            end
            [QN, UN, RN, TN, AN, WN] = deal(variables{:});  % Assign the modified values back to the original variables

            %%
            lqn = self.model.getStruct;
            Node = label(lqn.names);
            O = length(Node);
            NodeType = label(O,1);
            for o = 1:O
                switch lqn.type(o)
                    case LayeredNetworkElement.PROCESSOR
                        NodeType(o,1) = label({'Processor'});
                    case LayeredNetworkElement.TASK
                        if self.model.getStruct.isref(o)
                            NodeType(o,1) = label({'RefTask'});
                        else
                            NodeType(o,1) = label({'Task'});
                        end
                    case LayeredNetworkElement.ENTRY
                        NodeType(o,1) = label({'Entry'});
                    case LayeredNetworkElement.ACTIVITY
                        NodeType(o,1) = label({'Activity'});
                    case LayeredNetworkElement.CALL
                        NodeType(o,1) = label({'Call'});
                end
            end
            QLen = QN;
            QT = Table(Node,QLen);
            Util = UN;
            UT = Table(Node,Util);
            RespT = RN;
            RT = Table(Node,RespT);
            Tput = TN;
            TT = Table(Node,Tput);
            %SvcT = SN;
            %ST = Table(Node,SvcT);
            %ProcUtil = PN;
            %PT = Table(Node,ProcUtil);
            ResidT = WN;
            WT = Table(Node,ResidT);
            ArvR = AN;
            AT = Table(Node,ArvR);
            AvgTable = Table(Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput);%, ProcUtil, SvcT);
        end
    end

    methods % implemented in .m files
        runtime = runAnalyzer(self, options);
        [result, iterations] = parseXMLResults(self, filename);
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self);
        savedfname = plot(model);

        function allMethods = listValidMethods(self)
            %LISTVALIDMETHODS List valid solving methods for LQNS
            sn = self.model.getStruct();
            allMethods = {
                'default', 'lqns', 'srvn', 'exactmva', ...
                'srvn.exactmva', 'sim', 'lqsim', 'lqnsdefault'
                };
        end
    end

    methods (Static)


        function bool = isAvailable()
            %ISAVAILABLE Check if LQNS is available on the system path
            bool = true;
            if ispc
                [~, ret] = dos('lqns -V -H');
                if contains(ret, 'not recognized', 'IgnoreCase', true)
                    bool = false;
                    return;
                end
                if contains(ret, 'Version 5', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 4', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 3', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 2', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 1', 'IgnoreCase', true)
                    line_warning(mfilename, ...
                        'Unsupported LQNS version. LINE requires Version 6.0 or greater.');
                end
            else
                [~, ret] = unix('lqns -V -H');
                if contains(ret, 'command not found', 'IgnoreCase', true)
                    bool = false;
                    return;
                end
                if contains(ret, 'Version 5', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 4', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 3', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 2', 'IgnoreCase', true) || ...
                        contains(ret, 'Version 1', 'IgnoreCase', true)
                    line_warning(mfilename, ...
                        'Unsupported LQNS version. LINE requires Version 6.0 or greater.');
                end
            end
        end

        function [bool, featSupported] = supports(model)
            %SUPPORTS Check if the used features are supported
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFeatureSet;
            featSupported.setTrue({ ...
                'Sink', ...
                'Source', ...
                'Queue', ...
                'Coxian', ...
                'Erlang', ...
                'Exp', ...
                'HyperExp', ...
                'Buffer', ...
                'Server', ...
                'JobSink', ...
                'RandomSource', ...
                'ServiceTunnel', ...
                'SchedStrategy_PS', ...
                'SchedStrategy_FCFS', ...
                'ClosedClass' ...
                });

            bool = true;
            numLayers = model.getNumberOfLayers();
            for idx = 1:numLayers
                bool = bool && SolverFeatureSet.supports( ...
                    featSupported, featUsed{idx} ...
                    );
            end
        end

        function options = defaultOptions()
            %DEFAULTOPTIONS Return default options for SolverLQNS
            options = SolverOptions('LQNS');
        end

    end
end
