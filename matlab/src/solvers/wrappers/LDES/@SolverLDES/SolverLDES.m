classdef SolverLDES < NetworkSolver
    % SolverLDES LINE Discrete Event Simulator solver using SSJ library
    %
    % SolverLDES implements a discrete-event simulation solver that uses the SSJ
    % (Stochastic Simulation in Java) library to analyze queueing networks.
    % It supports open and closed networks with various service distributions,
    % scheduling strategies, and advanced node types.
    %
    % For LayeredNetwork (LQN) models, SolverLDES also supports LDES simulation.
    %
    % @brief Discrete-event simulation solver using SSJ library
    %
    % Example:
    % @code
    % solver = SolverLDES(model, 'samples', 1000000, 'seed', 23000);
    % solver.getAvg();  % Run LDES simulation
    % @endcode
    %
    % Passing an auxiliary solver as first optional argument warm-starts the
    % simulation from that solver's steady-state solution (see initFromSolver):
    % @code
    % solver = SolverLDES(model, SolverMVA(model), 'samples', 50000);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverLDES(model, varargin)
            % SOLVERLDES Create an LDES solver instance
            %
            % @brief Creates a Discrete Event Simulation solver
            % @param model Network model or path to a JMT file (.jsimg/.jsim/.jsimw/.jmva)
            % @param varargin Optional parameters (samples, seed, method, etc.)
            % @return self SolverLDES instance configured for simulation

            % Accept a JMT file path in place of a Network object
            if ischar(model) || isstring(model)
                model = JMT2LINE(char(model));
            end
            % An auxiliary solver passed as first optional argument requests a
            % warm start: its steady-state distribution decides the initial
            % simulation state (see initFromSolver).
            initSolver = [];
            if ~isempty(varargin) && isa(varargin{1}, 'NetworkSolver')
                initSolver = varargin{1};
                varargin(1) = [];
            end
            % LayeredNetwork (LQN) models are solved by the Java LDES backend
            % directly; attach the jline mirror before the superclass
            % constructor so the Network-specific initialization is skipped.
            if isa(model, 'LayeredNetwork') && isempty(model.obj)
                model.obj = JLINE.from_line_layered_network(model);
            end
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.options.lang = 'java';
            if isa(self.model, 'LayeredNetwork')
                % LayeredNetwork (LQN) models are solved by the Java LDES
                % ensemble backend; setLang builds self.obj (JLINE.SolverLDES).
                self.setLang();
            end
            % Regular Network models are fully JSON-mediated: the model is NOT
            % marshalled to a Java object (no self.obj, no model.obj), so
            % runAnalyzer uses the subprocess JSON path and repeated solvers keep
            % independent MATLAB-side transient handles (NetworkSolver.initHandles).
            if ~isempty(initSolver)
                self.initFromSolver(initSolver);
            end
        end

        function sn = getStruct(self)
            % QN = GETSTRUCT()

            % Get data structure summarizing the model
            if isa(self.model, 'LayeredNetwork')
                sn = self.model.getStruct();
            else
                sn = self.model.getStruct(true);
            end
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver

            allMethods = {'default'};
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSD>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % LDES is a discrete-event simulator; all methods are stochastic.
            bool = true;
        end
    end

    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink', 'Source', ...
                'Queue', 'Delay', ...
                'Fork', 'Join', 'Forker', 'Joiner', ... % Fork-Join node support
                'Place', 'Transition', ... % Petri net node support
                'QueueingPlace', ... % Queueing place (QPN embedded queue): FCFS/LCFS/SIRO/INF, renewal service
                'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage', ... % Petri net section support
                'Logger', 'LogTunnel', ... % Logger node support
                'Buffer', ...  % Finite buffer capacity support
                'Region', ...  % Finite capacity region support
                'Exp', 'Erlang', 'HyperExp', 'PH', 'APH', 'Coxian', 'Cox2', 'MAP', 'DMAP', 'MMAP', 'BMAP', 'MMPP2', 'ME', 'RAP', 'Immediate', 'Disabled', 'Replayer', 'Trace', ... % Trace is an alias of Replayer
                'Det', 'Uniform', 'Gamma', 'Pareto', 'Weibull', 'Lognormal', ... % Additional continuous distributions
                'Geometric', ... % Lattice-valued interarrival/service time on {1,2,...} (Geo/Geo/1 and slotted models)
                'Bernoulli', 'Binomial', 'Poisson', ... % Counting distributions; zero atom becomes an immediate interval (continuous mode only)
                'NHPP', ... % Piecewise-constant-intensity non-homogeneous Poisson process
                'Server', 'JobSink', 'RandomSource', ...
                'InfiniteServer', 'SharedServer', 'ServiceTunnel', 'DelayStation', ... % internal station-section markers
                'SchedStrategy_FCFS', 'SchedStrategy_INF', ...
                'SchedStrategy_HOL', ... % Priority scheduling (FCFS with priorities)
                'SchedStrategy_FCFSPRIO', ... % FCFS with priorities (non-preemptive)
                'SchedStrategy_PS', ... % Processor Sharing
                'SchedStrategy_DPS', ... % Discriminatory Processor Sharing
                'SchedStrategy_GPS', ... % Generalized Processor Sharing
                'SchedStrategy_LCFS', ... % Last Come First Served (non-preemptive)
                'SchedStrategy_LCFSPR', ... % LCFS Preemptive Resume
                'SchedStrategy_LCFSPI', ... % LCFS Preemptive Independent
                'SchedStrategy_FCFSPR', ... % FCFS Preemptive Resume
                'SchedStrategy_FCFSPI', ... % FCFS Preemptive Independent
                'SchedStrategy_LPS', ... % Longest Processing time first Shortest
                'SchedStrategy_SIRO', ...
                'SchedStrategy_SJF', 'SchedStrategy_LJF', ...
                'SchedStrategy_LEPT', ...
                'SchedStrategy_SEPT', ...
                'SchedStrategy_SRPT', ... % Shortest Remaining Processing Time (preemptive)
                'SchedStrategy_SRPTPRIO', ... % SRPT with priorities
                'SchedStrategy_PSJF', ... % Preemptive Shortest Job First
                'SchedStrategy_FB', ... % Feedback / Least Attained Service
                'SchedStrategy_LRPT', ... % Longest Remaining Processing Time
                'SchedStrategy_EXT', ...
                'SchedStrategy_POLLING', ... % Polling scheduling (GATED, EXHAUSTIVE, KLIMITED)
                'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO', 'SchedStrategy_GPSPRIO', ... % PS/DPS/GPS with priorities
                'SchedStrategy_LCFSPRIO', 'SchedStrategy_LCFSPRPRIO', 'SchedStrategy_LCFSPIPRIO', ... % LCFS priority variants
                'SchedStrategy_FCFSPRPRIO', 'SchedStrategy_FCFSPIPRIO', ... % FCFS preemptive priority variants
                'SchedStrategy_FSP', ... % Fair Sojourn Protocol (virtual PS finish time ranking)
                'SchedStrategy_PAS', ... % Pass-and-swap (order-independent) queue
                'SchedStrategy_OI', ... % Order-independent queue (PAS with empty swap graph)
                'SchedStrategy_EDD', 'SchedStrategy_EDF', 'SchedStrategy_SETF', ... % Deadline/elapsed-time disciplines
                'Router', 'Dispatcher', ... % Router node support (Dispatcher is the internal router section)
                'ClassSwitch', 'StatelessClassSwitcher', ... % Class switching node support
                'Cache', 'CacheClassSwitcher', ... % Cache node support with replacement policies (LRU, FIFO, Strict FIFO, RR)
                'RoutingStrategy_PROB', 'RoutingStrategy_RAND', ...
                'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN', ...
                'RoutingStrategy_JSQ', ... % Join the Shortest Queue
                'RoutingStrategy_SQ', ... % Power of K Choices routing
                'OpenClass', ...
                'ClosedClass', ...
                'SelfLoopingClass', ...
                'OpenSignal', ...           % G-network signal class in open networks
                'ClosedSignal', ...         % G-network signal class in closed networks
                'SignalType_NEGATIVE', ...
                'SignalType_REPLY', ...
                'SignalType_CATASTROPHE', ...
                'SignalBatchRemoval', ...   % Engine reads sn.signalremdist
                'SignalRemovalPolicy', ...  % Engine reads sn.signalrempolicy
                'LoadDependence', ... % Load-dependent service rates
                'ClassDependence', ... % Class-dependent service rate handles (setLimitedClassDependence)
                'JointDependence', ... % Joint-dependent (non-product-form) service rate handles (setJointDependence)
                'SetupDelayOff', ...  % Engine simulates the SETUP/DELAYOFF server states
                'Balking', ...        % Engine reads sn.balkingStrategy / balkingThresholds
                'Reneging', ...       % Engine collects renegingRate / avgRenegingWaitTime
                'Retrial', ...        % Engine collects retrialDropped and successful retries
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU',...
                'ReplacementStrategy_HLRU','ReplacementStrategy_CLIMB','ReplacementStrategy_QLRU'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            if isa(model, 'LayeredNetwork')
                % LayeredNetwork models are simulated by the Java LDES
                % backend (jline.solvers.ldes.SolverLDES), which validates
                % the LQN feature set at run time.
                bool = true;
                featSupported = SolverFeatureSet;
                return;
            end

            % Regular Network support
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverLDES.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('LDES');
        end

        function commonDir = getLdesCommonDir()
            % COMMONDIR = GETLDESCOMMONDIR()
            % Directory holding common/jline.jar and the optional native
            % GraalVM ldes binary. Resolved from the jline.jar entry on the
            % Java classpath (robust across checkouts and installs), falling
            % back to ascending from this class file to <root>/common.
            commonDir = '';
            try
                cp = javaclasspath('-all');
            catch
                cp = {};
            end
            for i = 1:numel(cp)
                [pdir, nm, ext] = fileparts(char(cp{i}));
                if strcmpi([nm ext], 'jline.jar')
                    commonDir = pdir;
                    return;
                end
            end
            % Fallback: <repo>/common relative to this file
            % (.../matlab/src/solvers/wrappers/LDES/@SolverLDES -> <repo> is 6 levels up).
            root = fileparts(mfilename('fullpath'));
            for k = 1:6
                root = fileparts(root);
            end
            cand = fullfile(root, 'common');
            if exist(cand, 'dir')
                commonDir = cand;
            end
        end

        function m = elfMachine(path)
            % M = ELFMACHINE(PATH)
            % ELF e_machine identifier of a binary (header offset 0x12, 2
            % bytes), or [] if PATH is not a readable ELF file.
            m = [];
            fid = fopen(path, 'r');
            if fid < 0
                return;
            end
            hdr = fread(fid, 20, '*uint8');
            fclose(fid);
            if numel(hdr) < 20
                return;
            end
            % Magic: 0x7F 'E' 'L' 'F'
            if ~(hdr(1) == 127 && hdr(2) == uint8('E') && hdr(3) == uint8('L') && hdr(4) == uint8('F'))
                return;
            end
            % EI_DATA (index 6, 1-based): 1 = little-endian, 2 = big-endian.
            % e_machine is a 2-byte field at offset 0x12 (indices 19,20).
            if hdr(6) == 1
                m = double(hdr(19)) + double(hdr(20)) * 256;
            else
                m = double(hdr(19)) * 256 + double(hdr(20));
            end
        end

        function m = hostElfMachine()
            % M = HOSTELFMACHINE()
            % Expected ELF e_machine for the current MATLAB host CPU, or [] if
            % unknown. MATLAB on Linux ships as glnxa64 (x86-64); glnxaa64
            % (aarch64) is mapped for completeness.
            switch computer('arch')
                case 'glnxa64'
                    m = 62;   % 0x3E EM_X86_64
                case 'glnxaa64'
                    m = 183;  % 0xB7 EM_AARCH64
                otherwise
                    m = [];
            end
        end

        function p = getLdesNativePath()
            % P = GETLDESNATIVEPATH()
            % Path to a runnable native LDES binary (common/ldes), or '' if
            % none is usable. Mirrors the Python-native selection: the native
            % binary is only used on Linux, and a binary whose ELF architecture
            % does not match the host CPU (e.g. an x86-64 build on an aarch64
            % host) is ignored so the caller can fall back to the in-process
            % JLINE (JVM) backend. If host or binary architecture cannot be
            % determined, the binary is used as a best effort.
            p = '';
            if ~(isunix && ~ismac)   % Linux only
                return;
            end
            commonDir = SolverLDES.getLdesCommonDir();
            if isempty(commonDir)
                return;
            end
            cand = fullfile(commonDir, 'ldes');
            if exist(cand, 'file') ~= 2
                return;
            end
            hostm = SolverLDES.hostElfMachine();
            binm = SolverLDES.elfMachine(cand);
            if ~isempty(hostm) && ~isempty(binm) && hostm ~= binm
                return;  % present but built for a different CPU architecture
            end
            p = cand;
        end

        function runners = getLdesRunners()
            % RUNNERS = GETLDESRUNNERS()
            % Ordered list of command prefixes that run the LDES engine on a
            % "solve ..." argument list, exchanging only JSON. The native GraalVM
            % binary (common/ldes) is tried first for fast startup; the full-JVM
            % "<java> -jar common/ldes.jar" is the fallback (same shaded engine).
            % The JVM fallback is needed because the AOT native binary lacks some
            % reflective/serialization features (e.g. the fork-join MMT transform
            % serializes jline.lang.Model). Returns a cellstr (possibly empty).
            runners = {};
            nativePath = SolverLDES.getLdesNativePath();
            if ~isempty(nativePath)
                runners{end+1} = sprintf('"%s"', nativePath);
            end
            commonDir = SolverLDES.getLdesCommonDir();
            if ~isempty(commonDir)
                ldesJar = fullfile(commonDir, 'ldes.jar');
                javaExe = SolverLDES.getJavaExe();
                if exist(ldesJar, 'file') == 2 && ~isempty(javaExe)
                    runners{end+1} = sprintf('"%s" -jar "%s"', javaExe, ldesJar);
                end
            end
        end

        function javaExe = getJavaExe()
            % JAVAEXE = GETJAVAEXE()
            % Resolve a Java launcher: LINE_JAVA, then JAVA_HOME/bin/java, then
            % the JRE bundled with MATLAB, then "java" on PATH. Returns '' if none
            % is found.
            javaExe = '';
            exeName = 'java';
            if ispc
                exeName = 'java.exe';
            end
            cands = {};
            envJava = getenv('LINE_JAVA');
            if ~isempty(envJava)
                cands{end+1} = envJava; %#ok<AGROW>
            end
            javaHome = getenv('JAVA_HOME');
            if ~isempty(javaHome)
                cands{end+1} = fullfile(javaHome, 'bin', exeName); %#ok<AGROW>
            end
            % JRE bundled with MATLAB (layout varies across releases).
            mlJre = fullfile(matlabroot, 'sys', 'java', 'jre', computer('arch'), 'jre', 'bin', exeName);
            cands{end+1} = mlJre; %#ok<AGROW>
            for i = 1:numel(cands)
                if exist(cands{i}, 'file') == 2
                    javaExe = cands{i};
                    return;
                end
            end
            % Last resort: rely on PATH resolution.
            [st, ~] = system(sprintf('%s -version', exeName));
            if st == 0
                javaExe = exeName;
            end
        end

    end
end
