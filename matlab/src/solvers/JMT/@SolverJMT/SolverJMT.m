classdef SolverJMT < NetworkSolver
    % A solver that interfaces the Java Modelling Tools (JMT) to LINE.
    %
    % Copyright (c) 2012-2023, Imperial College London
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
            % SELF = SOLVERJMT(MODEL, VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if ~Solver.isJavaAvailable
                line_error(mfilename,'SolverJMT requires the java command to be available on the system path.\n');
            end
            if ~SolverJMT.isAvailable
                line_warning(mfilename,'SolverJMT cannot located JMT.jar in the MATLAB path.\n');
            end
            self.simConfInt = 0.99;
            self.simMaxRelErr = 0.03;
            self.maxEvents = -1;
            self.setJMTJarPath(jmtGetPath);
        end
        
        [simDoc, section] = saveArrivalStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveBufferCapacity(self, simDoc, section, ind)
        [simDoc, section] = saveDropStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveGetStrategy(self, simDoc, section)
        [simDoc, section] = saveNumberOfServers(self, simDoc, section, ind)
        [simDoc, section] = savePreemptiveStrategy(self, simDoc, section, ind)
        [simDoc, section] = savePreemptiveWeights(self, simDoc, section, ind)
        [simDoc, section] = savePutStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveRoutingStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveServerVisits(self, simDoc, section)
        [simDoc, section] = saveServiceStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveClassSwitchStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveLogTunnel(self, simDoc, section, ind)
        [simDoc, section] = saveForkStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveJoinStrategy(self, simDoc, section, ind)
        [simElem, simDoc] = saveClasses(self, simElem, simDoc)
        [simElem, simDoc] = saveLinks(self, simElem, simDoc)
        [simElem, simDoc] = saveRegions(self, simElem, simDoc)
        [simElem, simDoc] = saveMetric(self, simElem, simDoc, handles)
        [simElem, simDoc] = saveMetrics(self, simElem, simDoc)
        [simElem, simDoc] = saveXMLHeader(self, logPath)
        
        fileName = getFileName(self)
        
        %Setter
        self = setJMTJarPath(self, path)
        
        % Getters
        out = getJMTJarPath(self)
        
        out = getFilePath(self)
        jsimwView(self, options)
        jsimgView(self, options)
        
        [outputFileName] = writeJSIM(self, sn, outputFileName)
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
        tranNodeStateAggr = sampleAggr(self, node, numSamples, markActivePassive);
        tranSysStateAggr = sampleSysAggr(self, numSamples, markActivePassive);
        
        %% Cdf methods
        [RD,log] = getCdfRespT(self, R);
        RD = getTranCdfRespT(self, R);
        RD = getTranCdfPassT(self, R);
    end
    
    methods (Static)
        
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','jsim','jmva','jmva.amva','jmva.mva','jmva.recal',...
                'jmva.comom','jmva.chow','jmva.bs','jmva.aql',...
                'jmva.lin','jmva.dmlin','jmva.ls'};
        end        

        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            
            bool = true;
            if isempty(which('JMT.jar'))
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
                'MAP',...
                'MMPP2',...
                'Normal',...
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
                'Linkage',...
                'Enabling', ...
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
                'SchedStrategy_LCFS',...
                'SchedStrategy_LCFSPR',...
                'SchedStrategy_SEPT',...
                'SchedStrategy_LEPT',...
                'SchedStrategy_SJF',...
                'SchedStrategy_LJF',...
                'RoutingStrategy_PROB',...
                'RoutingStrategy_RAND',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_WRROBIN',...
                'RoutingStrategy_KCHOICES',...
                'SchedStrategy_EXT',...
                'ClosedClass',...
                'OpenClass'});
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
            runtime = java.lang.Runtime.getRuntime();
            cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimg "',filename,'"'];
            system(cmd);
            %runtime.exec(cmd);
        end
        
        function jsimwOpen(filename)
            % JSIMWOPEN(FILENAME)
            
            runtime = java.lang.Runtime.getRuntime();
            cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',which(filename)];
            %system(cmd);
            runtime.exec(cmd);
        end
        
        dataSet = parseLogs(model, isNodeLogged, metric);
        [state, evtype, evclass, evjob] = parseTranState(fileArv, fileDep, nodePreload);
        [classResT, jobResT, jobResTArvTS, classResTJobID] = parseTranRespT(fileArv, fileDep);
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('JMT');
        end
        
        [outputFileName] = writeJMVA(sn, outputFileName, options)
    end
    
end
