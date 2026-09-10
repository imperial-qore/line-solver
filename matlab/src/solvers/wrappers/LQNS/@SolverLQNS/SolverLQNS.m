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
                    'Obtain them from their authors at: http://www.sce.carleton.ca/rads/lqns/\n' ...
                    'LINE ships no LQNS binary and does not redistribute one.\n\n' ...
                    'Alternatively, point LINE at a host that already runs LQNS:\n' ...
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

        function varargout = getAvgTable(self, varargin)
            % [AVGTABLE,QT,UT,RT,WT,TT] = GETAVGTABLE()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'avg', varargout{1});
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTable_impl(self)
            % GETAVGTABLE_IMPL Implementation of GETAVGTABLE; see the wrapper above.
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

        function bool = isStochasticMethod(self, method) %#ok<INUSL>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % The lqsim simulator is stochastic; the analytical lqns/srvn
            % methods are deterministic.
            bool = any(strcmpi(method, {'sim','lqsim'}));
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD)
            % The method-aware gate model.help and SolverAUTO ask; the rules
            % and their second caller, runAnalyzer, are in LQNS_METHOD_REFUSAL.
            [bool, reason] = lqns_method_refusal(self.model, method);
        end
    end

    methods (Static)


        function bool = hasLocalBinary()
            %HASLOCALBINARY True if the native lqns command is on the system path
            if ispc
                [~, ret] = dos('lqns -V -H');
                bool = ~contains(ret, 'not recognized', 'IgnoreCase', true);
            else
                [~, ret] = unix('lqns -V -H');
                bool = ~contains(ret, 'command not found', 'IgnoreCase', true);
            end
        end

        function bool = isAvailable()
            %ISAVAILABLE Check if a local lqns binary is installed
            %
            % LINE never runs LQNS from a container image: its licence is an
            % evaluation agreement that forbids redistribution, so the binary
            % must be one the user installed themselves. To exercise a
            % containerised LQNS in the test suite, put a shim on the PATH with
            % run-tests.sh --lqns-docker.
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

        function reasons = unsupportedLNConstructs(model)
            % REASONS = UNSUPPORTEDLNCONSTRUCTS(MODEL)
            % The LQN-level constructs neither lqns nor lqsim can model, named
            % one by one, or {} when the model carries none.
            %
            % THE FEATURE SET USED TO ADMIT THESE, and that is why the loss was
            % silent. `supports` below tests model.getUsedLangFeatures(), which
            % reports the features of the FLATTENED per-layer Networks and has no
            % vocabulary for a CacheTask, an ItemEntry or a task setup time; the
            % LQN-level construct therefore could not fail a check that never saw
            % it. The .lqnx writer now carries all three (see the LINE dialect in
            % writeXML.m), so lqns is handed a file it parses no further --
            % "Unexpected element <cache items=...>" -- which reports a
            % third-party syntax error for what is really a modelling limit of
            % the binary. Name the limit here instead. The writer stays correct
            % and unguarded: it is this feature set that was lying.
            reasons = {};
            for t = 1:numel(model.tasks)
                task = model.tasks{t};
                if isa(task, 'CacheTask')
                    reasons{end+1} = sprintf(['task ''%s'' is a CacheTask, and neither lqns nor ' ...
                        'lqsim models a cache'], task.name); %#ok<AGROW>
                end
                entries = task.entries;
                for e = 1:numel(entries)
                    if isa(entries(e), 'ItemEntry')
                        reasons{end+1} = sprintf(['entry ''%s'' on task ''%s'' is an ItemEntry, and ' ...
                            'neither lqns nor lqsim models an item reference stream'], ...
                            entries(e).name, task.name); %#ok<AGROW>
                    end
                end
                if isprop(task, 'setupTimeMean') && ~isempty(task.setupTimeMean) ...
                        && task.setupTimeMean > GlobalConstants.FineTol
                    reasons{end+1} = sprintf(['task ''%s'' declares a setup time (%g), which neither ' ...
                        'lqns nor lqsim charges'], task.name, task.setupTimeMean); %#ok<AGROW>
                end
                if isprop(task, 'delayOffTimeMean') && ~isempty(task.delayOffTimeMean) ...
                        && task.delayOffTimeMean > GlobalConstants.FineTol
                    reasons{end+1} = sprintf(['task ''%s'' declares a delay-off time (%g), which ' ...
                        'neither lqns nor lqsim charges'], task.name, task.delayOffTimeMean); %#ok<AGROW>
                end
            end
            % A queue-dependent service rate and a compatibility declaration are
            % both LQN-level and both invisible to the per-layer feature sets, for
            % the same reason the three above are. lqns has NO vocabulary for
            % either: its processors take a multiplicity and one of
            % {fcfs,hol,inf,ps,rand,pri}, its nine multiserver approximations are
            % all homogeneous, and the only class-differentiating mechanisms it
            % offers are priority and CFS group SHARES -- none of which is a
            % per-(server, class) eligibility. Answering for a homogeneous pool of
            % the same total size is a different system, so name the limit rather
            % than write a .lqnx that silently drops the structure.
            servers = [model.hosts(:); model.tasks(:)];
            for k = 1:numel(servers)
                elem = servers{k};
                if ~isa(elem, 'LayeredNetworkElement')
                    continue
                end
                if ismethod(elem, 'hasServerPools') && elem.hasServerPools()
                    reasons{end+1} = sprintf(['''%s'' declares heterogeneous server pools with a ' ...
                        'class-compatibility graph, which neither lqns nor lqsim models: their ' ...
                        'multiserver is homogeneous'], elem.name); %#ok<AGROW>
                end
                if isprop(elem, 'lldScaling') && ~isempty(elem.lldScaling)
                    reasons{end+1} = sprintf(['''%s'' declares a load-dependent service rate, which ' ...
                        'neither lqns nor lqsim models'], elem.name); %#ok<AGROW>
                end
                if isprop(elem, 'lcdScaling') && ~isempty(elem.lcdScaling)
                    reasons{end+1} = sprintf(['''%s'' declares a class-dependent service rate, which ' ...
                        'neither lqns nor lqsim models'], elem.name); %#ok<AGROW>
                end
                if isprop(elem, 'ljdScaling') && ~isempty(elem.ljdScaling)
                    reasons{end+1} = sprintf(['''%s'' declares a joint-dependent service rate, which ' ...
                        'neither lqns nor lqsim models'], elem.name); %#ok<AGROW>
                end
            end
        end

        function assertSupported(model, method)
            % ASSERTSUPPORTED(MODEL, METHOD)
            % Refuse a model this wrapper cannot serve, naming the construct.
            % Called before the .lqnx is written, so the answer is LINE's own and
            % not the binary's parse error. The rules are LQNS_METHOD_REFUSAL's,
            % the predicate the gate asks, so both speak one sentence.
            if nargin < 2
                method = '';
            end
            [ok, reason] = lqns_method_refusal(model, method);
            if ~ok
                line_error(mfilename, reason);
            end
        end

        function [bool, featSupported] = supports(model)
            %SUPPORTS Whether lqns can read this LayeredNetwork, whatever the method.
            % The method-neutral half of LQNS_METHOD_REFUSAL; see there for why
            % the per-layer feature comparison this used to make refused every
            % layered model. The second output keeps the signature the ensemble
            % callers expect: an LQN carries no flat feature set to compare.
            bool = lqns_method_refusal(model, '');
            featSupported = SolverFeatureSet;
        end

        function options = defaultOptions()
            %DEFAULTOPTIONS Return default options for SolverLQNS
            options = SolverOptions('LQNS');
        end

    end
end
