classdef JMTIO < handle
    % JMTIO Unified I/O handler for JMT model files
    %
    % Provides both read and write functionality for JMT XML format files
    % (JSIM and JMVA). Combines functionality previously split between
    % JMTXMLParser (write) and JMTResultParser (read).
    %
    % @brief Unified I/O for JMT simulation and analysis models
    %
    % Example:
    % @code
    % jmtio = JMTIO(model, options);
    % % Write
    % jmtio.writeJSIM(sn, outputFileName);
    % % Read
    % logData = JMTIO.parseLogs(model, isNodeLogged, metric);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        model           % Network model reference
        simConfInt      % Confidence interval for simulation
        simMaxRelErr    % Maximum relative error
        maxSamples      % Maximum samples
        maxEvents       % Maximum events
        maxSimulatedTime % Maximum simulated time
        seed            % Random seed
        fileName        % Output file name (without extension)
        options         % Solver options
        handles         % Metric handles from solver
    end

    properties (Constant)
        xsiNoNamespaceSchemaLocation = 'Archive.xsd';
        fileFormat = 'jsimg';
    end

    methods
        function self = JMTIO(model, options)
            % JMTIO Create a JMT I/O handler instance
            %
            % @brief Creates a unified I/O handler for JMT model files
            % @param model Network model to be converted to JMT format
            % @param options SolverOptions with simulation parameters
            % @return self JMTIO instance

            self.model = model;
            self.fileName = 'model';

            % Set defaults
            self.simConfInt = 0.99;
            self.simMaxRelErr = 0.03;
            self.maxEvents = -1;
            self.maxSamples = 10000;
            self.maxSimulatedTime = Inf;
            self.seed = 23000;

            % Override with options if provided
            if nargin >= 2 && ~isempty(options)
                self.options = options;
                if isfield(options, 'samples') || isprop(options, 'samples')
                    self.maxSamples = options.samples;
                end
                if isfield(options, 'seed') || isprop(options, 'seed')
                    self.seed = options.seed;
                end
                % Parse confidence interval from options
                [confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);
                if confintEnabled
                    self.simConfInt = confintLevel;
                end
            else
                self.options = struct('seed', self.seed);
            end
        end

        function sn = getStruct(self)
            % GETSTRUCT Get the network structure
            %
            % @return sn Network structure from model
            sn = self.model.getStruct(true);
        end

        % External method declarations (in separate files)
        fileName = getFileName(self)
        [exportClasses, cacheClasses] = getExportableClasses(self)
        [simElem, simDoc] = saveXMLHeader(self, logPath)
        [simElem, simDoc] = saveClasses(self, simElem, simDoc)
        [simDoc, section] = saveArrivalStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveBufferCapacity(self, simDoc, section, ind)
        [simDoc, section] = saveCacheStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveClassSwitchStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveDelayOffStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveDropRule(self, simDoc, section, ind)
        [simDoc, section] = saveDropStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveEnablingConditions(self, simDoc, section, ind)
        [simDoc, section] = saveFCRMetrics(self, simElem, simDoc)
        [simDoc, section] = saveFiringOutcomes(self, simDoc, section, ind)
        [simDoc, section] = saveFiringPriorities(self, simDoc, section, ind)
        [simDoc, section] = saveFiringWeights(self, simDoc, section, ind)
        [simDoc, section] = saveForkStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveGetStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveImpatience(self, simDoc, section, ind)
        [simDoc, section] = saveInhibitingConditions(self, simDoc, section, ind)
        [simDoc, section] = saveJoinStrategy(self, simDoc, section, ind)
        [simElem, simDoc] = saveLinks(self, simElem, simDoc)
        [simDoc, section] = saveLogTunnel(self, simDoc, section, ind)
        [simElem, simDoc] = saveMetric(self, simElem, simDoc, handles)
        [simElem, simDoc] = saveMetrics(self, simElem, simDoc)
        [simDoc, section] = saveModeNames(self, simDoc, section, ind)
        [simDoc, section] = saveNumberOfServers(self, simDoc, section, ind)
        [simDoc, section] = saveNumbersOfServers(self, simDoc, section, ind)
        [simDoc, section] = savePlaceCapacities(self, simDoc, section, ind)
        [simDoc, section] = savePreemptiveStrategy(self, simDoc, section, ind)
        [simDoc, section] = savePreemptiveWeights(self, simDoc, section, ind)
        [simDoc, section] = savePutStrategies(self, simDoc, section, ind)
        [simDoc, section] = savePutStrategy(self, simDoc, section, ind)
        [simElem, simDoc] = saveRegions(self, simElem, simDoc)
        [simDoc, section] = saveRoutingStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveServerVisits(self, simDoc, section)
        [simDoc, section] = saveServiceStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveSwitchoverStrategy(self, simDoc, section, ind)
        [simDoc, section] = saveTimingStrategies(self, simDoc, section, ind)
        [simDoc, section] = saveTotalCapacity(self, simDoc, section, ind)

        outputFileName = writeJSIM(self, sn, outputFileName)
    end

    methods (Static)
        outputFileName = writeJMVA(sn, outputFileName, options)

        % Parse JMT log files for transient metrics
        logData = parseLogs(model, isNodeLogged, metric)

        % Parse transient state data from arrival/departure logs
        [state, evtype, evclass, evjob] = parseTranState(fileArv, fileDep, nodePreload)

        % Parse transient response time data from arrival/departure logs
        [classResT, jobResT, jobResTArvTS, classResTJobID] = parseTranRespT(fileArv, fileDep)
    end
end
