classdef SolverCTMC < NetworkSolver
    % A solver based on continuous-time Markov chain (CTMC) formalism.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverCTMC(model,varargin)
            % SELF = SOLVERCTMC(MODEL,VARARGIN)            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
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
            
        [infGen, eventFilt, synchInfo, stateSpace, nodeStateSpace] = getSymbolicGenerator(self, invertSymbol, primeNumbers)
        [infGen, eventFilt, synchInfo] = getInfGen(self, options)        
        [infGen, eventFilt, synchInfo] = getGenerator(self, options)
        
        tstate = sampleSys(self, numevents)
        sampleAggr = sampleAggr(self, node, numSamples)
        
        function MCTMC = getMarkedCTMC(self, options)        
            % MCTMC = GETMARKEDCTMC(options)

            if nargin < 2
                [infGen, eventFilt, synchInfo] = self.getInfGen();    
            else
                [infGen, eventFilt, synchInfo] = getInfGen(self, options);    
            end
            
            MCTMC = MarkedCTMC(infGen, eventFilt, synchInfo);
        end

        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
    end
           
    methods (Static)

        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','gpu'};
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Source','Sink',...
                'ClassSwitch','DelayStation','Queue',...
                'MAP','APH','MMPP2','PH','Coxian','Erlang','Exp','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Cache','CacheClassSwitcher', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_GPS',...
                'SchedStrategy_SIRO','SchedStrategy_SEPT',...
                'SchedStrategy_LEPT','SchedStrategy_FCFS',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'SchedStrategy_LCFSPR',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','OpenClass','Replayer'});
        end
        
        function [bool, featSupported, featUsed] = supports(model)
            % [BOOL, FEATSUPPORTED, FEATUSED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverCTMC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('CTMC');
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
    end
end
