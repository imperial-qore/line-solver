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
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
            if ~Solver.isJavaAvailable
                line_error(mfilename,'SolverJMT requires the java command to be available on the system path.\n');
            end
            if ~SolverJMT.isAvailable
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
            self.setJMTJarPath(jmtGetPath);
            % Initialize XML parser for model serialization
            self.xmlParser = JMTIO(model, self.options);
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
        
        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default','jsim','replication','jmva','jmva.amva','jmva.mva','jmva.recal',...
                'jmva.comom','jmva.chow','jmva.bs','jmva.aql',...
                'jmva.lin','jmva.dmlin'};
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
                'SchedStrategy_SRPT',...
                'SchedStrategy_LEPT',...
                'SchedStrategy_SJF',...
                'SchedStrategy_LJF',...
                'SchedStrategy_LPS',...
                'SchedStrategy_POLLING',...
                'RoutingStrategy_PROB',...
                'RoutingStrategy_RAND',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_WRROBIN',...
                'RoutingStrategy_KCHOICES',...
                'SchedStrategy_EXT',...
                'ClosedClass','SelfLoopingClass',...
                'OpenClass',...
                'Cache', 'CacheClassSwitcher', ...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', ...
                'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU'});
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
