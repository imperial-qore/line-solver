classdef SolverJMT < NetworkSolver
    % SolverJMT Java Modelling Tools interface for simulation analysis
    %
    % SolverJMT provides an interface to the Java Modelling Tools (JMT) suite
    % for discrete event simulation of queueing networks. It handles model
    % translation, simulation execution, and results parsing for both JMVA
    % (analytical) and JSIM (simulation) engines within JMT.
    %
    % @brief Interface to Java Modelling Tools for simulation and analysis
    %
    % Example:
    % @code
    % solver = SolverJMT(model, 'samples', 10000, 'seed', 1);
    % solver.getAvg();  % Run simulation via JMT
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    %Private properties
    properties %(GetAccess = 'private', SetAccess='private')
        jmtPath;
        filePath;
        fileName;
        maxSimulatedTime;
        maxSamples;
        maxEvents;
        seed;
        simConfInt;
        simMaxRelErr;
        xmlParser;  % JMTIO instance for XML generation
    end
    
    %Constants
    properties (Constant)
        xsiNoNamespaceSchemaLocation = 'Archive.xsd';
        fileFormat = 'jsimg';
        jsimgPath = '';
    end
    
    % PUBLIC METHODS
    methods
        
        %Constructor
        function self = SolverJMT(model, varargin)
            % SOLVERJMT Create a JMT solver instance
            %
            % @brief Creates a Java Modelling Tools solver for the given model
            % @param model Network model to be analyzed via JMT
            % @param varargin Optional parameters (samples, seed, maxSimulatedTime, etc.)
            % @return self SolverJMT instance configured for JMT simulation

            % An auxiliary solver passed as first optional argument requests a
            % warm start: its steady-state solution decides the station
            % preload of the simulation (see NetworkSolver.initFromSolver).
            initSolver = [];
            if ~isempty(varargin) && isa(varargin{1}, 'NetworkSolver')
                initSolver = varargin{1};
                varargin(1) = [];
            end
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
            % A missing JVM is no longer fatal: jmtRun can still reach JMT
            % through options.rest_url or the Docker image, and it raises an
            % actionable error if neither is usable. The jar is resolved only
            % on the local path, so a JVM-free host does not download 50MB it
            % will never load.
            hasJava = Solver.isJavaAvailable;
            if ~hasJava
                line_warning(mfilename,'The java command is not on the system path. SolverJMT will use a JMT REST server (options.rest_url) or the JMT Docker image instead.\n');
            elseif ~SolverJMT.isAvailable
                line_warning(mfilename,'SolverJMT cannot locate JMT.jar in the common folder.\n');
            end
            % Set confidence interval from options (default: 0.99 for 99% confidence)
            [confintEnabled, confintLevel] = Solver.parseConfInt(self.options.confint);
            if confintEnabled
                self.simConfInt = confintLevel;
            else
                self.simConfInt = 0.99; % default when disabled
            end
            self.simMaxRelErr = 0.03;
            self.maxEvents = -1;
            if hasJava
                self.setJMTJarPath(jmtGetPath);
            end
            % Initialize XML parser for model serialization
            self.xmlParser = JMTIO(model, self.options);
            if ~isempty(initSolver)
                self.initFromSolver(initSolver);
            end
        end

        % XML generation delegated to JMTIO (see @JMTIO/)
        % Use self.xmlParser.writeJSIM() for JSIM model generation

        fileName = getFileName(self)
        
        %Setter
        self = setJMTJarPath(self, path)
        
        % Getters
        out = getJMTJarPath(self)
        
        out = getFilePath(self)
        jsimwView(self, options)
        jsimgView(self, options)
        view(self, options)

        [outputFileName] = writeJSIM(self, sn, outputFileName)

        % Backend dispatch: local JVM, JMT REST server, or Docker image
        [status, cmdout] = jmtRun(self, mode, fname, seed, options)
        [status, cmdout] = jmtSolveRest(self, restUrl, mode, fname, seed, options)

        [result, parsed] = getResults(self)
        [result, parsed] = getResultsJSIM(self)
        [result, parsed] = getResultsJMVA(self)
        
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
    end
    
    %Private methods.
    methods (Access = 'private')
        out = getJSIMTempPath(self)
        out = getJMVATempPath(self)
    end
    
    %Private methods.
    methods (Access = 'protected')
        bool = hasAvgResults(self)
    end
    
    
    methods (Access = 'public')
        getProbNormConstAggr(self); % jmva
        %% StateAggr methods
        Pr = getProbAggr(self, node, state_a);
        [Pi_t, SSnode_a] = getTranProbAggr(self, node);
        probSysStateAggr = getProbSysAggr(self);
        tranNodeStateAggr = sampleAggr(self, node, numEvents, markActivePassive);
        tranSysStateAggr = sampleSysAggr(self, numEvents, markActivePassive);
        
        %% Cdf methods
        [RD,log] = getCdfRespT(self, R);
        RD = getTranCdfRespT(self, R);
        RD = getTranCdfPassT(self, R);
        
        function bool = supportsTransientAnalysis(self) %#ok<MANU>
            % Transient averages are available (simulation restricted to options.timespan).
            bool = true;
        end

        function method = resolveMethod(self, options) %#ok<INUSL>
            % METHOD = RESOLVEMETHOD(OPTIONS)
            % 'default' runs the JSIM steady-state arm, or the 'replication'
            % transient arm on a finite timespan: the same switch runAnalyzer
            % makes, stated here so runAnalyzerChecks gates 'default' as the
            % arm that will actually run and names it when it refuses.
            method = options.method;
            if strcmpi(method, 'default')
                if isfield(options,'timespan') && numel(options.timespan) >= 2 ...
                        && isfinite(options.timespan(2))
                    method = 'replication';
                else
                    method = 'jsim';
                end
            end
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default','jsim','replication','jmva','jmva.amva','jmva.mva','jmva.recal',...
                'jmva.comom','jmva.chow','jmva.bs','jmva.aql',...
                'jmva.lin','jmva.dmlin'};
        end

        function reason = unsupportedMethodReason(self, method) %#ok<INUSL>
            % REASON = UNSUPPORTEDMETHODREASON(METHOD)
            %
            % The forwarding address for 'jmva.ls', JMVA's Logistic Sampling
            % (SolverAlgorithm.MONTE_CARLO_LOGISTIC, algType name 'Logistic
            % Sampling'). The engine still implements it; LINE withdrew the
            % name on 2025-06-04 (b345d7d4e), which dropped it from
            % listValidMethods and commented out the writeJMVA arm that emits
            % the algType. Without this the name gate answers with the flat
            % "unsupported by this solver" and the caller cannot tell a
            % withdrawn method from a typo, nor find the LINE-native estimator
            % of the same normalizing constant. Asks nothing of the model, so
            % runAnalyzerChecks can call it before the struct is built.
            reason = '';
            if ischar(method) && any(strcmp(method, {'jmva.ls','jmt.jmva.ls'}))
                reason = sprintf(['The ''%s'' method (JMVA Logistic Sampling) is withdrawn: ' ...
                    'writeJMVA no longer emits the ''Logistic Sampling'' algType, so the name ' ...
                    'would silently run exact MVA. Use SolverNC(model,''method'',''ls'') for ' ...
                    'the logistic-sampling normalizing constant, or ''jmva.recal''/''jmva.comom'' ' ...
                    'for an exact one.'], method);
            end
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSL>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % Simulation-based methods (default, jsim, replication) return
            % stochastic estimates. The analytical JMVA methods do not,
            % except for its sampling-based variants (e.g. jmva.ls).
            methodNames = regexp(lower(method), '[./]', 'split');
            if any(strcmp(methodNames, 'jmva'))
                bool = any(ismember(methodNames, {'ls','mci','imci','sampling'}));
            else
                bool = true;
            end
        end
        function featSupported = getMethodFeatureSet(self, method)
            % SolverJMT drives TWO ENGINES, and they accept different models.
            %
            % 'default', 'jsim' and 'replication' run the JSIM SIMULATOR, whose
            % envelope is getFeatureSet() below: caches, fork-join, Petri nets,
            % finite capacity regions, impatience, heterogeneous server pools,
            % every discipline the writer emits. The 'jmva.*' names run the
            % JMVA ANALYTICAL engine, which reads a document carrying only a
            % station type (delay / load-independent / load-dependent), a
            % per-chain service demand, a per-chain visit count, the class
            % populations or arrival rates and a reference station. Everything
            % else in the model is DROPPED by writeJMVA, so declaring the JSIM
            % envelope for jmva was a promise the writer could not keep: on a
            % three-class LRU cache model all eight closed-form jmva methods
            % returned an entirely zero table with no error, jmva.mva labelled
            % 'exact' among them.
            %
            % Defining this is also what lets NetworkSolver.supportsModelMethod
            % name the offending features: with no method feature set it falls
            % back to the coarse supports(model), which returns an empty reason.
            %
            % A non-Network model (e.g. a LayeredNetwork) has no
            % getUsedLangFeatures, so it keeps the coarse path and any
            % structural checks or redirects that operate on such models.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            if strncmpi(method, 'jmva', 4)
                featSupported = SolverJMT.getJMVAFeatureSet();
                if jmtJmvaIsClosedOnly(method)
                    % RECAL, CoMoM, Chow, Bard-Schweitzer, AQL, Linearizer and
                    % De Souza-Muntz Linearizer are closed-network algorithms:
                    % JMT answers an open or a mixed model with "The selected
                    % solver cannot handle open classes" and a load-dependent
                    % one with the matching refusal. Exact MVA, which 'jmva'
                    % and 'jmva.mva' select, serves both.
                    featSupported.setFalse({'OpenClass','LoadDependence'});
                    % the eight closed-form algorithms are single-server ones
                    % (jmtMethodRefusal words it); exact MVA carries the count
                    featSupported.setFalse({'MultiServer'});
                end
                return;
            end
            featSupported = SolverJMT.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Structural gate for what no registry name can state: the finite
            % timespan the 'replication' arm integrates over, the single-server
            % restriction of the closed-form JMVA algorithms (a server count is
            % not a declared feature), immediate feedback, the mean-only JMVA
            % document (see jmtMethodRefusal), and the one feature JMT admits
            % in a RE-ENCODED form only. Limited load dependence has no
            % representation of its own in either JMT document:
            % saveNumberOfServers turns it into a server count and writeJMVA
            % into the matching <ldstation>, so alpha(n) = min(n,c) with an
            % integer c is written exactly and any other scaling would be
            % solved at a service rate JMT never saw.
            if isa(self.model, 'Network')
                sn = self.model.getStruct();
                % The same predicate writeJMVA and the replication arm ask, so
                % this gate and those runs cannot answer differently.
                structural = jmtMethodRefusal(sn, method, self.getOptions());
                if ~isempty(structural)
                    bool = false;
                    reason = structural;
                    return;
                end
                if sn_has_load_dependence(sn)
                    for ist=1:min(sn.nstations, size(sn.lldscaling,1))
                        alpha = sn.lldscaling(ist,:);
                        c = max(alpha);
                        if all(alpha == 1)
                            continue
                        end
                        if c ~= round(c) || c < 1 || ...
                                any(abs(alpha - min(1:numel(alpha), c)) > GlobalConstants.Zero)
                            bool = false;
                            reason = sprintf(['Station %d uses a load-dependent scaling that is ' ...
                                'not the multiserver encoding alpha(n) = min(n,c): JMT has no ' ...
                                'representation for it, since both the JSIM and the JMVA writer ' ...
                                'carry the scaling as a server count, and the model would be ' ...
                                'solved at the nominal service rate. Use SolverCTMC, SolverNC, ' ...
                                'SolverMVA or SolverSSA, which read sn.lldscaling directly.'], ist);
                            return;
                        end
                    end
                end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
        end

    end
    
    methods (Static)

        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            
            bool = true;
            try
                jmt_path = jmtGetPath();
                jmt_jar = fullfile(jmt_path, 'JMT.jar');
                if ~exist(jmt_jar, 'file')
                    bool = false;
                end
            catch
                bool = false;
            end
        end
        
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
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
                'JoinPartial',... % quorum join, written out as a jmt PartialJoin
                ... A variable forking level: saveForkStrategy turns
                ... isSimplifiedFork off and writes the per-branch entries, so
                ... jmt reads the counts, the probabilities and the degree
                ... distribution rather than sending one job down every link.
                'ForkFanoutVector',...
                'ForkFanoutRandom',...
                'ForkBranchProbability',...
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
                'MAP',...
                'MMPP2',...
                ... % 'Normal' is NOT declared: MATLAB has no ProcessType.NORMAL, so a
                ... % Normal cannot even be attached as a service or arrival process
                ... % ("Unrecognized process type"), and neither saveServiceStrategy
                ... % nor saveArrivalStrategy has an arm for it. Dead registry surface
                ... % for services here until both exist.
                'PH',...
                'Pareto',...
                'Weibull',...                                
                'Replayer',...
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
                'Buffer', ...
                ... % Finite capacity regions: JMT is the reference FCR engine --
                ... % every other solver's FCR refusal message points users here.
                'Region', ...
                'Linkage',...
                'Enabling', ...
                'Inhibiting', ...
                'Timing', ...
                'Firing', ...
                'Storage', ...
                'Place', ...
                'Transition', ...
                'SchedStrategy_INF',...
                'SchedStrategy_PS',...
                'SchedStrategy_DPS',...
                'SchedStrategy_FCFS',...
                'SchedStrategy_GPS',...
                'SchedStrategy_SIRO',...
                'SchedStrategy_HOL',...
                'SchedStrategy_PSPRIO',...
                'SchedStrategy_DPSPRIO',...
                'SchedStrategy_GPSPRIO',...
                'SchedStrategy_LCFS',...
                'SchedStrategy_LCFSPR',...
                'SchedStrategy_LCFSPRIO',...
                'SchedStrategy_LCFSPRPRIO',...
                ...% The preemptive family savePutStrategy.m already emits
                ...% (lines 63-91) and writeJSIM.m already promotes to
                ...% PreemptiveServer. Declared here because the writer supports
                ...% them: withholding the names refused a model this wrapper
                ...% exports correctly, which is how the C++ port came to declare
                ...% six scheduling names the reference did not.
                'SchedStrategy_LCFSPI',...
                'SchedStrategy_LCFSPIPRIO',...
                'SchedStrategy_FCFSPR',...
                'SchedStrategy_FCFSPI',...
                'SchedStrategy_FCFSPRPRIO',...
                'SchedStrategy_FCFSPIPRIO',...
                ...% Earliest Due Date and Earliest Deadline First are REAL in
                ...% JMT 1.2.x, not placeholders: EDDStrategy.class and
                ...% EDFStrategy.class ship in JMT.jar and order the buffer by
                ...% the due dates saveClassSoftDeadlines writes as
                ...% <classSoftDeadlines>. Verified 2026-09-05 by running
                ...% jmt.commandline.Jmt on a two-class M/M/1 of equal rates:
                ...% with deadlines 5.0 and 1.0 the tight class got R = 3.69
                ...% against 5.56, which FCFS cannot produce. A class WITHOUT a
                ...% due date aborts the run, so jmtDeadlineRefusal gates it.
                'SchedStrategy_EDD',...
                'SchedStrategy_EDF',...
                'SchedStrategy_SEPT',...
                'SchedStrategy_SRPT',...
                'SchedStrategy_SRPTPRIO',...
                'SchedStrategy_LEPT',...
                'SchedStrategy_SJF',...
                'SchedStrategy_LJF',...
                'SchedStrategy_LPS',...
                'SchedStrategy_POLLING',...
                'RoutingStrategy_PROB',...
                'RoutingStrategy_RAND',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_WRROBIN',...
                'RoutingStrategy_JSQ',...
                'RoutingStrategy_SQ',...
                'SchedStrategy_EXT',...
                'ClosedClass','SelfLoopingClass',...
                'OpenClass',...
                'Cache', 'CacheClassSwitcher', ...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', ...
                'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU', ...
                ... % Limited load dependence reaches JMT only as a SERVER COUNT:
                ... % saveNumberOfServers exports max(nservers,max(alpha)) and
                ... % writeJMVA writes the matching <ldstation>. That is exact for
                ... % alpha(n) = min(n,c) and for nothing else, so
                ... % supportsModelMethod refuses any other scaling by name.
                'LoadDependence', ...
                'SetupDelayOff', ... % exported as delayOffTime/setUpTime (saveDelayOffStrategy)
                'ServerParallelism', ... % exported as classParallelism (saveClassParallelism)
                ... % Heterogeneous server pools: the type names, the servers per
                ... % type and the compatibility matrix are exported as
                ... % serverTypesNames / serverTypesNumOfServers /
                ... % serverTypesCompatibilities (saveServerTypeNames,
                ... % saveServersPerType, saveServerCompatibilities), so jsim
                ... % simulates the pools rather than a station of the same total
                ... % size.
                'HeteroServers', ...
                'Reneging', ...      % exported as an Impatience/Reneging strategy (saveImpatience)
                'Balking', ...       % exported as an Impatience/Balking strategy (saveImpatience)
                ... % Queue.setRetrial: writeJSIM picks the retrial Queue constructor and
                ... % saveRetrialDistributions writes the per-class orbit delay. Batch
                ... % arrivals (Source.setArrivalBatch) are NOT declared: saveArrivalStrategy
                ... % has no batch element, so the stream would be written single.
                'Retrial', ...
                ... % c-server stations (saveNumberOfServers) and finite buffers
                ... % (saveBufferCapacity with the drop rule): the JSIM writer
                ... % exports both, and jmtStationCapRefusal keeps refusing the
                ... % buffers JMT would answer unconstrained (closed WAITQ, BBS,
                ... % RSRD, RETRIAL_WITH_LIMIT); the JMVA set withdraws the buffer.
                'MultiServer', 'FiniteCapacity'});
        end

        function featSupported = getJMVAFeatureSet()
            % FEATSUPPORTED = GETJMVAFEATURESET()
            %
            % What the JMVA ANALYTICAL engine accepts, which is much less than
            % the JSIM simulator above. The envelope is derived from the writer
            % rather than guessed: writeJMVA emits, per station, a
            % <delaystation>, a <listation> or an <ldstation>, a per-chain
            % <servicetime> and a per-chain <visit>, and at model level the
            % closed populations, the open arrival rates and the reference
            % station. NOTHING ELSE IN THE MODEL REACHES JMVA, so a construct
            % whose whole effect is not carried by (station type, demand,
            % visits, population) would be solved away silently.
            %
            % The names dropped from the JSIM set, and why each is dropped:
            %   Cache and the replacement strategies -- no cache element
            %     exists; the hit/miss split and the cache state are lost. This
            %     is the case that returned an all-zero table.
            %   Fork/Join and the fan-out names -- no fork element exists, and
            %     a visit ratio cannot express the join synchronization.
            %   Place/Transition and the Petri-net sections -- no counterpart.
            %   Region -- JMVA has no finite capacity region.
            %   Reneging/Balking -- no impatience element; the abandonment
            %     would simply not happen.
            %   SetupDelayOff, ServerParallelism, HeteroServers -- each is a
            %     server-side attribute the JMVA document has no slot for.
            %   the non-BCMP disciplines -- the writer emits NO discipline at
            %     all, so a priority, weighted, size-based or limited-sharing
            %     station would be solved as an ordinary load-independent one.
            %     Only the four BCMP station types survive the encoding, which
            %     is the same line SolverNC and SolverMVA draw.
            %   the state-dependent routings (RROBIN, WRROBIN, JSQ, SQ) -- the
            %     document carries mean visit counts, which is not what makes a
            %     join-the-shortest-queue model behave as it does.
            % The DISTRIBUTIONS are deliberately kept: JMVA consumes a mean
            % service demand, so any renewal law with a finite mean is
            % admissible, exactly as it is for SolverMVA and SolverNC.
            featSupported = SolverJMT.getFeatureSet();
            featSupported.setFalse({...
                'Cache','CacheClassSwitcher', ...
                'ReplacementStrategy_RR','ReplacementStrategy_FIFO', ...
                'ReplacementStrategy_SFIFO','ReplacementStrategy_LRU', ...
                'Fork','Join','Forker','Joiner','JoinPartial', ...
                'ForkFanoutVector','ForkFanoutRandom','ForkBranchProbability', ...
                'Place','Transition','Enabling','Inhibiting','Timing', ...
                'Firing','Storage', ...
                'Region', ...
                'Reneging','Balking','Retrial', ...
                'SetupDelayOff','ServerParallelism','HeteroServers', ...
                'SchedStrategy_DPS','SchedStrategy_GPS','SchedStrategy_HOL', ...
                'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO', ...
                'SchedStrategy_LCFSPI','SchedStrategy_LCFSPIPRIO', ...
                'SchedStrategy_LCFSPRIO','SchedStrategy_LCFSPRPRIO', ...
                'SchedStrategy_FCFSPR','SchedStrategy_FCFSPI', ...
                'SchedStrategy_FCFSPRPRIO','SchedStrategy_FCFSPIPRIO', ...
                'SchedStrategy_EDD','SchedStrategy_EDF', ...
                'SchedStrategy_SEPT','SchedStrategy_LEPT', ...
                'SchedStrategy_SJF','SchedStrategy_LJF', ...
                'SchedStrategy_SRPT','SchedStrategy_SRPTPRIO', ...
                'SchedStrategy_LPS','SchedStrategy_POLLING', ...
                'RoutingStrategy_RROBIN','RoutingStrategy_WRROBIN', ...
                'RoutingStrategy_JSQ','RoutingStrategy_SQ', ...
                ... % the JMVA document has no capacity element at all
                'FiniteCapacity'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverJMT.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function jsimgOpen(filename)
            % JSIMGOPEN(FILENAME)
            
            [path] = fileparts(filename);
            if isempty(path)
                filename=[pwd,filesep,filename];
            end
            cmd = [line_java_cmd(mfilename),' -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimg "',filename,'"'];
            system(cmd);
        end

        function jsimwOpen(filename)
            % JSIMWOPEN(FILENAME)

            cmd = [line_java_cmd(mfilename),' -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',which(filename),'"'];
            system(cmd);
        end
        
        % Parse methods delegated to JMTResultParser (see @JMTResultParser/)
        % For backward compatibility, wrapper methods exist in @SolverJMT/
        dataSet = parseLogs(model, isNodeLogged, metric);
        [state, evtype, evclass, evjob] = parseTranState(fileArv, fileDep, nodePreload);
        [classResT, jobResT, jobResTArvTS, classResTJobID] = parseTranRespT(fileArv, fileDep);

        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('JMT');
        end
        
        [outputFileName] = writeJMVA(sn, outputFileName, options)
    end
    
end
