classdef SolverLQNS < NetworkSolver
    % A solver that interfaces the LQNS toolset to LINE.
    %
    % A LayeredNetwork is solved by lqns or lqsim. A flat Network is solved by
    % qnsolver, the product-form MVA solver of the same distribution, through
    % the 'qns' methods; this is the arrangement SolverJMT has for JMVA.
    %
    % Methods on a LayeredNetwork:
    %   default, lqns, srvn, exactmva, srvn.exactmva, sim, lqsim, lqns.default
    % Methods on a Network, each naming the multiserver approximation:
    %   qns (and default) - rolia, the qnsolver default here
    %   qns.conway - Conway (1989), extending the multinomial all-servers-busy probability of de Souza e Silva and Muntz (Perform. Eval. 7(3), 1987)
    %   qns.rolia - Rolia (PhD thesis, Toronto, 1992) as used in the method of layers (Rolia and Sevcik, IEEE TSE 21(8), 1995), in the per-class Rolia-Franks form of Franks (PhD thesis, Carleton, 1999)
    %   qns.zhou - arrival-theorem binomial (AB) approximation, S. Zhou (M.A.Sc. thesis, Carleton, 2021) and Zhou and Woodside (ICPE Companion 2022)
    %   qns.suri - Suri, Sahu and Vernon (IERC 2007)
    %   qns.reiser - Reiser and Lavenberg (J. ACM 27(2), 1980) load-dependent MVA, see also Reiser (Perform. Eval. 1, 1981)
    %   qns.schmidt - Schmidt (Perform. Eval. 29(4), 1997)
    % A product-form or open Network goes to qnsolver; a closed non-product-form
    % one is converted by QN2LQN and goes to lqns (runAnalyzerNetwork).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverLQNS(model, varargin)
            % SELF = SOLVERLQNS(MODEL, VARARGIN)
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            % A Network is served by qnsolver, and by lqns only on the closed
            % non-product-form branch, so the lqns binary is not a precondition
            % of constructing one: that branch reports its own absence.
            if isa(model, 'LayeredNetwork') && ~SolverLQNS.isAvailable() && ~self.options.config.remote
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
            if isa(self.model, 'LayeredNetwork')
                sn = self.model.getStruct();
            else
                sn = self.model.getStruct(false); % doesn't need initial state
            end
        end

        function [QN,UN,RN,TN,AN,WN] = getAvgLayered(self)
            %GETAVGLAYERED Per-element metrics of a LayeredNetwork, from lqns/lqsim
            % NetworkSolver.getAvg, which is sealed, takes this branch for a
            % LayeredNetwork; a Network goes through runAnalyzerNetwork.
            [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self);
        end

        function [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTableLayered(self)
            % GETAVGTABLELAYERED Per-element table of a LayeredNetwork.
            % NetworkSolver.getAvgTable, which is sealed and records the table
            % it returns (see LineResultRecorder), takes this branch for one.
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
        [runtime, analyzer] = runAnalyzer(self, options);
        [runtime, analyzer] = runAnalyzerNetwork(self, options);
        [result, iterations] = parseXMLResults(self, filename);
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self);
        savedfname = plot(model);

        function allMethods = listValidMethods(self)
            %LISTVALIDMETHODS List valid solving methods for LQNS
            % The 'qns' names reach qnsolver and take a flat Network; the others
            % reach lqns/lqsim and take a LayeredNetwork.
            if isa(self.model, 'LayeredNetwork')
                allMethods = {
                    'default', 'lqns', 'srvn', 'exactmva', ...
                    'srvn.exactmva', 'sim', 'lqsim', 'lqns.default'
                    };
            else
                allMethods = SolverLQNS.qnsMethods();
            end
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSL>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % The lqsim simulator is stochastic; the analytical lqns/srvn
            % methods are deterministic.
            bool = any(strcmpi(method, {'sim','lqsim'}));
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % All qns methods share one feature envelope, getFeatureSet.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model), which returns an empty reason, so
            % the gate could only report "features not supported" without
            % saying which ones. A LayeredNetwork has no getUsedLangFeatures and
            % keeps the coarse path.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverLQNS.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD)
            % On a LayeredNetwork, the method-aware gate model.help and
            % SolverAUTO ask; the rules and their second caller, runAnalyzer,
            % are in LQNS_METHOD_REFUSAL. On a Network, the qnsolver gate below.
            if isa(self.model, 'LayeredNetwork')
                if SolverLQNS.isQnsMethod(method)
                    bool = false;
                    reason = sprintf(['SolverLQNS: method ''%s'' runs qnsolver, which solves a flat ' ...
                        'Network; a LayeredNetwork takes lqns or lqsim (%s).'], char(method), ...
                        strjoin(self.listValidMethods(), ', '));
                    return
                end
                [bool, reason] = lqns_method_refusal(self.model, method);
                return
            end
            % Structural gate, on top of the feature envelope, for what no
            % registry name can state.
            %
            % A BINDING FINITE BUFFER first. NOTHING under the qns methods reads
            % sn.cap or sn.classcap -- neither the JMVA document qnsolver reads
            % nor the LQN that QN2LQN writes has a buffer -- so a capped station
            % was solved as an unbounded one and the table reported the
            % unconstrained answer under this solver's name. SolverMVA,
            % SolverNC, SolverAG and SolverFLD gate the same way through the
            % same helper. Without it SolverAUTO.listValidMethods offered all
            % eight 'qns' method names on the BAS-blocking model of cqn_bas_blocking.
            %
            % THEN THE PATH. runAnalyzerNetwork sends a product-form or an open model
            % to qnsolver through writeJMVA and a closed non-product-form one
            % to lqns through QN2LQN, and the two carry different things.
            % On the qnsolver path runAnalyzerNetwork maps the method name onto
            % options.config.multiserver one to one, and 'qnsolver -m' knows
            % conway, reiser, rolia and zhou only, so 'suri' and 'schmidt' die
            % on a multiserver station: QNS_MULTISERVER_REFUSAL is the predicate
            % SOLVER_QNS asks at run time, so it is asked here too. (The former
            % comment here claimed the config value stayed 'default' under
            % those names and Conway answered; runAnalyzerNetwork sets it from the
            % name before SOLVER_QNS runs, so it does not.) The JMVA document
            % also has no fork or join and carries mean demands only, which
            % JMTMETHODREFUSAL states for the engine, exactly as writeJMVA asks
            % it before writing. On the lqns path QN2LQN writes AND forks and
            % joins as activity precedences and passes the full service law, so
            % the same model is served there and no path rule applies.
            %
            % IMMEDIATE FEEDBACK is refused on BOTH paths, in this solver's own
            % words (QNS_IMMFEED_REFUSAL, also asked by runAnalyzerNetwork): neither
            % document can keep a self-looping job on its server.
            % STRUCTURAL PREDICATE FIRST, as SolverBA does and for the same
            % reason: it names the station, the cap and the way out, where the
            % feature envelope can only say "(feature: FiniteCapacity)". Once
            % FiniteCapacity became a registry name on 2026-09-05 the base gate
            % started answering first and the useful sentence became unreachable.
            if isa(self.model, 'Network')
                [bool, reason] = NetworkSolver.checkBindingCapacity(self.model, 'SolverLQNS');
                if ~bool
                    return
                end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if ~bool || ~isa(self.model, 'Network')
                return
            end
            sn = self.model.getStruct();
            reason = qns_immfeed_refusal(sn);
            if ~isempty(reason)
                bool = false;
                return
            end
            if self.model.hasProductFormSolution() || self.model.hasOpenClasses()
                [bool, reason] = qns_multiserver_refusal(sn, SolverLQNS.qnsMultiserver(method));
                if ~bool
                    return
                end
                reason = jmtMethodRefusal(sn, method, self.getOptions(), 'jmva');
                bool = isempty(reason);
            end
        end
    end

    methods (Static)

        function names = qnsMethods()
            % NAMES = QNSMETHODS()
            % The method names on a flat Network: 'default' and 'qns' take the
            % rolia multiserver approximation, 'qns.<name>' names another one.
            names = {'default', 'qns', 'qns.conway', 'qns.rolia', 'qns.zhou', ...
                'qns.suri', 'qns.reiser', 'qns.schmidt'};
        end

        function bool = isQnsMethod(method)
            % BOOL = ISQNSMETHOD(METHOD) True for 'qns' and every 'qns.<name>'.
            m = lower(char(method));
            bool = strcmp(m, 'qns') || strncmp(m, 'qns.', 4);
        end

        function ms = qnsMultiserver(method)
            % MS = QNSMULTISERVER(METHOD)
            % The multiserver approximation a Network method selects: the part
            % after 'qns.', or 'default' for 'default' and 'qns'.
            m = lower(char(method));
            if strncmp(m, 'qns.', 4)
                ms = m(5:end);
            else
                ms = 'default';
            end
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            %
            % What BOTH paths of runAnalyzerNetwork carry, i.e. the envelope of the
            % 'qns' methods on a flat Network (a LayeredNetwork has none). qnsolver reads the JMVA
            % document writeJMVA emits (a station type, a mean demand and a
            % visit count per chain) and lqns reads the LQN QN2LQN writes (a
            % host per station with its multiplicity and discipline, a task,
            % an entry per class, the service law, OR/AND forks from the
            % routing), so a construct is declared only when neither drops it.
            %
            % WITHDRAWN, and why: the Petri-net names, which writeJMVA ignores
            % and QN2LQN has no arm for, so a Place or Transition vanished; MAP
            % and MMPP2, whose correlation the JMVA document flattens to a rate
            % (undeclared, needsMapEnv now solves them through the random
            % environment image instead); every discipline outside INF, PS and
            % FCFS, since the JMVA document names no discipline and the .lqnx
            % schema knows no 'dps', 'siro', 'sept' or 'lcfs' (lqns rejects the
            % file), while HOL kept its name but lost the class priorities on
            % both paths; and the state-dependent routings (RROBIN, WRROBIN,
            % SQ), which mean visit counts cannot express. Fork/Join stay: the
            % lqns path writes them as precedences, and the qnsolver path,
            % which cannot, is refused by supportsModelMethod. The distributions
            % stay too, by the JMVA policy: a mean is admissible, and the mean-
            % only refusal at a FCFS station is jmtMethodRefusal's to state.
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink',...
                'Source',...
                'Router',...
                'ClassSwitch',...
                'Delay',...
                'DelayStation',...
                'Queue',...
                'Fork',...
                'Join',...
                'Forker',...
                'Joiner',...
                'Logger',...
                'Coxian',...
                'Cox2',...
                'APH',...
                'Erlang',...
                'Exp',...
                'HyperExp',...
                'Det',...
                'Gamma',...
                'Lognormal',...
                'Normal',...
                'PH',...
                'Pareto',...
                'Weibull',...
                'Replayer',...
                'Trace',...
                'Uniform',...
                'StatelessClassSwitcher',...
                'InfiniteServer',...
                'SharedServer',...
                'Buffer',...
                'Dispatcher',...
                'Server',...
                'JobSink',...
                'RandomSource',...
                'ServiceTunnel',...
                'LogTunnel',...
                'SchedStrategy_INF',...
                'SchedStrategy_PS',...
                'SchedStrategy_FCFS',...
                'RoutingStrategy_PROB',...
                'RoutingStrategy_RAND',...
                'SchedStrategy_EXT',...
                'ClosedClass',...
                'OpenClass',...
                ... % c-server stations: the JMVA document carries the count as an
                ... % <ldstation> and the LQN as a host multiplicity; 'suri' and
                ... % 'schmidt' refuse one on the qnsolver path, which stays
                ... % structural (QNS_MULTISERVER_REFUSAL). FiniteCapacity is NOT
                ... % declared: neither document has a buffer (supportsModelMethod).
                'MultiServer'});
        end

        function bool = hasQnsolver()
            %HASQNSOLVER True if the native qnsolver command is on the system path
            if ispc
                [~, ret] = dos('qnsolver -H');
                bool = ~contains(ret, 'not recognized', 'IgnoreCase', true);
            else
                [~, ret] = unix('qnsolver -H');
                bool = ~contains(ret, 'command not found', 'IgnoreCase', true);
            end
        end

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
            %SUPPORTS Whether this solver can read MODEL, whatever the method.
            % A LayeredNetwork: the method-neutral half of LQNS_METHOD_REFUSAL;
            % see there for why the per-layer feature comparison this used to
            % make refused every layered model. The second output keeps the
            % signature the ensemble callers expect: an LQN carries no flat
            % feature set to compare. A Network: the qns feature envelope.
            if isa(model, 'LayeredNetwork')
                bool = lqns_method_refusal(model, '');
                featSupported = SolverFeatureSet;
                return
            end
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverLQNS.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function options = defaultOptions()
            %DEFAULTOPTIONS Return default options for SolverLQNS
            options = SolverOptions('LQNS');
        end

    end
end
