classdef SolverCTMC < NetworkSolver
    % Continuous-Time Markov Chain solver for exact state-space analysis
    %
    % Implements exact analysis of queueing networks via CTMC formulation.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverCTMC(model,varargin)
            % SOLVERCTMC Create a CTMC solver instance
            %
            % @brief Creates a CTMC solver for exact Markov chain analysis
            % @param model Network model to be analyzed via CTMC formulation
            % @param varargin Optional parameters (cutoff, method, etc.)
            % @return self SolverCTMC instance configured for exact analysis

            % Redirect a LayeredNetwork to the solvers that accept one. This
            % must precede the NetworkSolver constructor: that constructor calls
            % model.getAvgHandles(), which a LayeredNetwork does not implement,
            % so the redirect that used to live in runAnalyzer was unreachable
            % and the user saw "Unrecognized method 'getAvgHandles'" instead.
            if isa(model, 'LayeredNetwork')
                line_error(mfilename, ['This model is a LayeredNetwork. Use SolverLN ', ...
                    '(iterative) or SolverLQNS (external LQNS) instead of SolverCTMC.\n']);
            end
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
        end
        
        runtime = runAnalyzer(self, options)
        Pnir = getProb(self, node, state)
        Pn = getProbSys(self)
        Pnir = getProbAggr(self, ist)
        
        Pn = getProbSysAggr(self)
        [Pi_t, SSsysa] = getTranProbSysAggr(self)   
        [Pi_t, SSnode_a] = getTranProbAggr(self, node)
        [Pi_t, SSsys] = getTranProbSys(self)     
        [Pi_t, SSnode] = getTranProb(self, node)
        RD = getCdfRespT(self, R)
        RD = getCdfSysRespT(self)
        
        [stateSpace,nodeStateSpace] = getStateSpace(self, options)
        stateSpaceAggr = getStateSpaceAggr(self)

        % Reward computation methods
        [Rt, t, names] = getTranReward(self, rewardName)
        [R, names] = getAvgReward(self)
        [V, t, names, stateSpace] = runRewardAnalyzer(self)
        
        function [state_space, local_states] = stateSpace(self)
            % STATESPACE Kotlin-style alias for getStateSpace
            if nargout <= 1
                state_space = self.getStateSpace();
            else
                [state_space, local_states] = self.getStateSpace();
            end
        end
        
        function Q = generator(self)
            % GENERATOR Kotlin-style alias for getGenerator
            Q = self.getGenerator();
        end
            
        [infGen, eventFilt, synchInfo, stateSpace, nodeStateSpace] = getSymbolicGenerator(self, invertSymbol)
        [infGen, eventFilt, synchInfo] = getInfGen(self, options)        
        [infGen, eventFilt, synchInfo] = getGenerator(self, options)
        
        tstate = sampleSys(self, numevents)
        sampleAggr = sampleAggr(self, node, numEvents)
        
        function MCTMC = getMarkedCTMC(self, options)        
            % MCTMC = GETMARKEDCTMC(options)

            if nargin < 2
                [infGen, eventFilt, synchInfo] = self.getInfGen();    
            else
                [infGen, eventFilt, synchInfo] = getInfGen(self, options);    
            end
            
            MCTMC = MarkedMarkovProcess(infGen, eventFilt, synchInfo);
        end

        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            % QRF (Quadratic/Linear Reduction) LP-based bounds were moved to
            % SolverBA (methods qr/lr/qrf.*); SolverCTMC no longer serves them.
            allMethods = {'default','gpu'};
        end
        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % All CTMC methods share the solver-level feature envelope.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model), which returns an empty reason, so
            % the gate could only report "features not supported" without
            % saying which ones.
            %
            % A non-Network model (e.g. a LayeredNetwork) has no
            % getUsedLangFeatures, so it keeps the coarse path and any
            % structural checks or redirects that operate on such models.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverCTMC.getFeatureSet();
        end


    end
           
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Source','Sink',...
                'ClassSwitch','Delay','DelayStation','Queue','Router',...
                'MAP','APH','MMPP2','MMAP','PH','Coxian','Erlang','Exp','HyperExp','ME',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                ... % Finite capacity regions: CTMC represents them exactly (the
                ... % blocked-job overflow buffer is part of the chain), verified by
                ... % convergence to the simulators as the cutoff grows.
                'Region', ...
                'Cache','CacheClassSwitcher', ...
                'CacheRetrieval', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_GPS',...
                'SchedStrategy_SIRO','SchedStrategy_SEPT',...
                'SchedStrategy_LEPT','SchedStrategy_FCFS',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'SchedStrategy_LCFSPR','SchedStrategy_LCFSPRPRIO','SchedStrategy_FCFSPRPRIO',...
                'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO',...
                'SchedStrategy_LPS',...
                'SchedStrategy_PAS',...
                'SchedStrategy_OI',...
                'SchedStrategy_POLLING',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_WRROBIN',...
                'RoutingStrategy_JSQ',...
                'RoutingStrategy_SQ',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO','ReplacementStrategy_LRU',...
                'ReplacementStrategy_HLRU','ReplacementStrategy_CLIMB','ReplacementStrategy_QLRU',...
                'ClosedClass','SelfLoopingClass','OpenClass','Replayer',...
                'OpenSignal','ClosedSignal',...
                'SignalType_NEGATIVE','SignalType_CATASTROPHE','SignalType_REPLY',...
                'SignalBatchRemoval','SignalRemovalPolicy',...
                'Place','Transition','Linkage','Enabling','Inhibiting','Timing','Firing','Storage',...
                'Fork','Join','Forker','Joiner',...
                'Balking','Reneging','Retrial','Breakdown',...
                'LoadDependence','ClassDependence','JointDependence'});
        end
        
        function [bool, featSupported, featUsed] = supports(model)
            % [BOOL, FEATSUPPORTED, FEATUSED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverCTMC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('CTMC');
        end
        
        function printInfGen(Q,SS)
            % PRINTINFGEN(Q,SS)
            
            SS=full(SS);
            Q=full(Q);
            for s=1:size(SS,1)
                for sp=1:size(SS,1)
                    if Q(s,sp)>0
                        line_printf('\n%s->%s: %f',mat2str(SS(s,:)),mat2str(SS(sp,:)),double(Q(s,sp)));
                    end
                end
            end
            line_printf('\n');
        end
        
        function printEventFilt(sync,D,SS,myevents)
            % PRINTEVENTFILT(SYNC,D,SS,MYEVENTS)
            
            if nargin<4 %~exist('events','var')
                myevents = 1:length(sync);
            end
            SS=full(SS);
            for e=myevents
                D{e}=full(D{e});
                for s=1:size(SS,1)
                    for sp=1:size(SS,1)
                        if D{e}(s,sp)>0
                            line_printf('\n%s-- %d: (%d,%d) => (%d,%d) -->%s: %f',mat2str(SS(s,:)),e,sync{e}.active{1}.node,sync{e}.active{1}.class,sync{e}.passive{1}.node,sync{e}.passive{1}.class,mat2str(SS(sp,:)),double(D{e}(s,sp)));
                        end
                    end
                end
            end
        end

        function libs = getLibrariesUsed(sn, options) %#ok<INUSD>
            % GETLIBRARIESUSED Get list of external libraries used by CTMC solver
            libs = {};
        end

        function backend = symbolicBackend(self)
            % BACKEND = SYMBOLICBACKEND(SELF)
            %
            % Value of options.config.symbolic, or 'auto' when the solver
            % carries no such option (an older options struct, or none at all).
            backend = 'auto';
            try
                if isprop(self, 'options') && isfield(self.options, 'config') && ...
                        isfield(self.options.config, 'symbolic')
                    backend = self.options.config.symbolic;
                end
            catch
            end
        end

        function entries = symbolicEntries(infGenTerms, symbols, invertSymbol, n)
            % ENTRIES = SYMBOLICENTRIES(INFGENTERMS, SYMBOLS, INVERTSYMBOL, N)
            %
            % Renders the per-event coefficient matrices as expression strings,
            % e.g. '2*x1 - 3*x2'. This is the wire form the symbolic backend
            % consumes and the exact format the JAR's getSymbolicEntry emits,
            % which is what lets the two be compared at all.
            entries = cell(n, n);
            for i = 1:n
                for j = 1:n
                    s = '';
                    for e = 1:numel(symbols)
                        if isempty(symbols{e})
                            continue
                        end
                        c = infGenTerms{e}(i, j);
                        if c == 0
                            continue
                        end
                        if isempty(s)
                            if c < 0
                                s = '-';
                            end
                        else
                            if c < 0
                                s = [s, ' - ']; %#ok<AGROW>
                            else
                                s = [s, ' + ']; %#ok<AGROW>
                            end
                        end
                        absC = abs(c);
                        if invertSymbol
                            s = [s, SolverCTMC.symbolicCoeff(absC), '/', symbols{e}]; %#ok<AGROW>
                        else
                            if absC ~= 1
                                s = [s, SolverCTMC.symbolicCoeff(absC), '*']; %#ok<AGROW>
                            end
                            s = [s, symbols{e}]; %#ok<AGROW>
                        end
                    end
                    if isempty(s)
                        s = '0';
                    end
                    entries{i, j} = s;
                end
            end
        end

        function s = symbolicCoeff(c)
            % S = SYMBOLICCOEFF(C)  Coefficient as text the backend reads exactly.
            if c == round(c) && abs(c) < 1e15
                s = sprintf('%d', round(c));
            else
                s = sprintf('%.17g', c);
            end
        end
    end
end
