function runtime = runAnalyzer(self, options)
% TSIM = RUN()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart=tic;

if nargin<2 %%~exist('options','var')
    options = self.getOptions;
end

line_debug('JMT analyzer starting: lang=%s, samples=%d, seed=%d', options.lang, options.samples, options.seed);

self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

switch options.lang
    case 'java'
        line_debug('Using Java JMT implementation, calling JLINE SolverJMT');
        jmodel = LINE2JLINE(self.model);
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverJMT(jmodel, options);
        T0=tic;
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);
        runtime = toc(T0);
        CN = [];
        XN = [];
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        WN = reshape(WN',R,M)';
        AN = reshape(AN',R,M)';
        TN = reshape(TN',R,M)';
        lG = NaN;
        lastiter = NaN;
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        self.result.solverSpecific.sn = JLINE.from_jline_struct(jmodel);
        self.result.Prob.logNormConstAggr = lG;
        return
    case 'matlab'
        line_debug('Using MATLAB JMT implementation');

        if ~isfield(options,'verbose')
            options.verbose = 0;
        end

        if ~isfield(options,'force')
            options.force = false;
        end

        if ~isfield(options,'keep')
            options.keep = false;
        end

        if self.enableChecks && ~self.supports(self.model)
            %    if options.verbose
            line_error(mfilename,'This model contains features not supported by the solver.');
            %    end
            %    runtime = toc(T0);
            %    return
        end

        if ~isfield(options,'samples')
            options.samples = 1e4; % default: this is the samples / measure, not the total number of simulation events, which can be much larger.
        elseif options.samples < 5e3
            %if ~strcmpi(options.method,'jmva.ls')
                line_warning(mfilename,'JMT requires at least 5000 samples for each metric, the current value is %d. Starting the simulation with 5000 samples.\n', options.samples);
            %end
            options.samples = 5e3;
        end

        if ~isfield(options,'verbose')
            options.verbose = 0;
        end

        if ~isfield(options,'keep')
            options.verbose = false;
        end

        if ~isfield(options,'seed')
            options.seed = randi([1,1e6]);
        end
        self.seed = options.seed;

        if ~isfield(options,'timespan')
            options.timespan = [0,Inf];
        else
            self.maxSimulatedTime = options.timespan(2);
        end

        if ~self.model.hasInitState
            self.model.initDefault;
        end

        self.maxSamples = options.samples;
        % Sync xmlParser maxSamples with solver (ensures corrected sample count is used in XML)
        self.xmlParser.maxSamples = options.samples;
        sn = self.getStruct;

        % Auto-detect transient mode: switch to 'replication' if finite timespan with default method
        if strcmpi(options.method, 'default') && isfield(options, 'timespan') && isfinite(options.timespan(2))
            options.method = 'replication';
            line_debug('Finite timespan detected, automatically using replication method for transient simulation');
        end

        switch options.method
            case {'jsim','default'}
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using JSIM discrete-event simulation\n');
                end
                line_debug('Using JSIM method (discrete-event simulation)');
                fname = self.writeJSIM(sn);
                cmd = ['java -cp "',getJMTJarPath(self),filesep,'JMT.jar" jmt.commandline.Jmt sim "',fname,'" -seed ',num2str(options.seed),' --illegal-access=permit'];
                if options.verbose
                    line_printf('JMT model: %s\n',fname);
                end
                if options.verbose == VerboseLevel.DEBUG
                    line_printf('JMT command: %s\n',cmd);
                end
                [status, cmdoutput] = system(cmd);
                runtime = toc(Tstart);
                if status ~= 0 && options.verbose
                    line_printf('\nJMT command failed with status %d. Output:\n%s\n', status, cmdoutput);
                end
                self.getResults;
                if ~options.keep && isfolder(getFilePath(self))
                    rmdir(getFilePath(self),'s');
                end
                %if options.verbose
                %    line_printf('\nJMT analysis (seed: %d) completed. Runtime: %f seconds.\n',options.seed,runtime);
                %end
                self.setAvgResults(self.result.Avg.Q,self.result.Avg.U,self.result.Avg.R,self.result.Avg.T,self.result.Avg.A,self.result.Avg.W,[],[],runtime,options.method,1);

                % Set cache hit/miss ratios from JMT Cache Hit Rate metric
                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        hitclass = sn.nodeparam{ind}.hitclass;
                        nclasses = length(hitclass);
                        actualhitprob = zeros(1, nclasses);
                        actualmissprob = zeros(1, nclasses);

                        % Check if JMT returned Cache Hit Rate metrics
                        fieldName = sprintf('node%d', ind);
                        if isfield(self.result, 'CacheHitRate') && isfield(self.result.CacheHitRate, fieldName)
                            actualhitprob = self.result.CacheHitRate.(fieldName);
                            actualmissprob = 1 - actualhitprob;
                        end

                        % Set the results on the Cache node
                        self.model.nodes{ind}.setResultHitProb(actualhitprob);
                        self.model.nodes{ind}.setResultMissProb(actualmissprob);

                        % Also update the cached sn struct directly
                        sn.nodeparam{ind}.actualhitprob = actualhitprob;
                        sn.nodeparam{ind}.actualmissprob = actualmissprob;
                        self.model.sn.nodeparam{ind}.actualhitprob = actualhitprob;
                        self.model.sn.nodeparam{ind}.actualmissprob = actualmissprob;
                    end
                end
            case {'replication'}
                line_debug('Using replication method (transient simulation)');
                options = self.getOptions;
                initSeed = self.options.seed;
                initTimeSpan = self.options.timespan;
                self.options.timespan(1) = self.options.timespan(2);
                if isfield(options,'timespan') && isfinite(options.timespan(2))
                    tu = [];
                    validReplications = 0;
                    for it=1:options.iter_max
                        self.options.seed = initSeed + it -1;
                        TranSysStateAggr{it} = sampleSysAggr(self);
                        % Skip replications with empty or invalid time vectors
                        if isempty(TranSysStateAggr{it}.t) || ~isvector(TranSysStateAggr{it}.t)
                            line_warning(mfilename, 'Replication %d produced empty/invalid time series, skipping.', it);
                            continue;
                        end
                        validReplications = validReplications + 1;
                        if isempty(tu)
                            tu = TranSysStateAggr{it}.t;
                        else
                            % we need to limit the time series at the minimum
                            % as otherwise the predictor of the state cannot
                            % take into account constraints that exist on the
                            % state space
                            tumax = min(max(tu),max(TranSysStateAggr{it}.t));
                            tu = union(tu, TranSysStateAggr{it}.t);
                            tu = tu(tu<=tumax);
                        end
                    end
                    if validReplications == 0
                        line_error(mfilename, 'No valid replications produced. Cannot compute transient averages.');
                        return;
                    end
                    QNt = cellzeros(sn.nstations, sn.nclasses, length(tu), 2);
                    UNt = cellzeros(sn.nstations, sn.nclasses, length(tu), 2);
                    TNt = cellzeros(sn.nstations, sn.nclasses, length(tu), 2);
                    M = sn.nstations;
                    K = sn.nclasses;
                    for jst=1:M
                        for r=1:K
                            QNt{jst,r}(:,2) = tu;
                            UNt{jst,r}(:,2) = tu;
                            TNt{jst,r}(:,2) = tu;
                            for it=1:options.iter_max
                                % Skip invalid replications
                                if isempty(TranSysStateAggr{it}.t) || ~isvector(TranSysStateAggr{it}.t)
                                    continue;
                                end
                                qlenAt_t = interp1(TranSysStateAggr{it}.t, TranSysStateAggr{it}.state{jst}(:,r), tu,'previous');
                                avgQlenAt_t = qlenAt_t;
                                %avgQlenAt_t = cumsum(qlenAt_t .*[0;diff(tu)])./tu;
                                avgQlenAt_t(isnan(avgQlenAt_t))=0;
                                QNt{jst,r}(:,1) = QNt{jst,r}(:,1) + (1/validReplications) * avgQlenAt_t;
                            end
                            for it=1:options.iter_max
                                % Skip invalid replications
                                if isempty(TranSysStateAggr{it}.t) || ~isvector(TranSysStateAggr{it}.t)
                                    continue;
                                end
                                if isfinite(sn.nservers(jst))
                                    occupancyAt_t = interp1(TranSysStateAggr{it}.t, min(TranSysStateAggr{it}.state{jst}(:,r),sn.nservers(jst)), tu,'previous')/sn.nservers(jst);
                                else % if delay we use queue-length
                                    occupancyAt_t = interp1(TranSysStateAggr{it}.t, TranSysStateAggr{it}.state{jst}(:,r), tu,'previous');
                                end
                                avgOccupancyAt_t = occupancyAt_t;
                                %avgOccupancyAt_t = cumsum(occupancyAt_t .*[0;diff(tu)])./tu;
                                avgOccupancyAt_t(isnan(avgOccupancyAt_t))=0;
                                UNt{jst,r}(:,1) = UNt{jst,r}(:,1) + (1/validReplications) * avgOccupancyAt_t;
                            end
                            %                         for it=1:options.iter_max
                            %                             departures = [0;diff(TranSysStateAggr{it}.state{j}(:,r))];
                            %                             departures(departures>0) = 0;
                            %                             departuresAt_t = abs(interp1(TranSysStateAggr{it}.t, cumsum(departures), tu, 'previous'));
                            %                             avgDeparturesAt_t = departuresAt_t./tu;
                            %                             avgDeparturesAt_t(isnan(avgDeparturesAt_t))=0;
                            %                             TNt{j,r}(:,1) = TNt{j,r}(:,1) + (1/options.iter_max) * avgDeparturesAt_t;
                            %                         end
                            if isfinite(sn.nservers(jst))
                                TNt{jst,r}(:,1) = UNt{jst,r}(:,1) * sn.nservers(jst) * sn.rates(jst,r);
                            else
                                TNt{jst,r}(:,1) = UNt{jst,r}(:,1) * sn.rates(jst,r);
                            end
                        end
                    end
                    runtime = toc(Tstart);
                    RNt = [];
                    CNt = [];
                    XNt = [];
                    self.setTranAvgResults(QNt,UNt,RNt,TNt,CNt,XNt,runtime);
                    self.result.Tran.Avg.U = UNt;
                    self.result.Tran.Avg.T = TNt;
                    self.result.Tran.Avg.Q = QNt;
                end
                self.options.seed = initSeed;
                self.options.timespan = initTimeSpan;
                self.result.('solver') = getName(self);
                self.result.runtime = runtime;
                %if options.verbose
                %   line_printf('\nJMT analysis (seed: %d) completed. Runtime: %f seconds.\n',options.seed,runtime);
                %end
            case {'jmva','jmva.amva','jmva.mva','jmva.recal','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmva.ls',...
                    'jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin','jmt.jmva.ls'}
                line_debug('Using JMVA method: %s', options.method);
                fname = self.writeJMVA(sn, getJMVATempPath(self), self.options);
                cmd = ['java -cp "',getJMTJarPath(self),filesep,'JMT.jar" jmt.commandline.Jmt mva "',fname,'" -seed ',num2str(options.seed),' --illegal-access=permit'];
                if options.verbose
                    line_printf('JMT model: %s\n',fname);
                end
                if options.verbose == VerboseLevel.DEBUG
                    line_printf('JMT command: %s\n',cmd);
                end
                [status, cmdoutput] = system(cmd);
                runtime = toc(Tstart);
                if status ~= 0 && options.verbose
                    line_printf('\nJMT command failed with status %d. Output:\n%s\n', status, cmdoutput);
                end
                self.getResults;
                if ~options.keep && isfolder(getFilePath(self))
                    rmdir(getFilePath(self),'s');
                end
                %if options.verbose
                %    line_printf('\nJMT analysis (method: %d) completed. Runtime: %f seconds.\n',options.method,runtime);
                %end
                 sn = self.model.getStruct();

                 % Compute average arrival rate at steady-state
                 TH = getAvgTputHandles(self);
                 AN = sn_get_arvr_from_tput(sn, self.result.Avg.T, TH);
                 self.setAvgResults(self.result.Avg.Q,self.result.Avg.U,self.result.Avg.R,self.result.Avg.T,AN,self.result.Avg.W,[],[],runtime,options.method,1);
            otherwise
                line_warning(mfilename,'This solver does not support the specified method. Setting to default.\n');
                self.options.method  = 'default';
                runAnalyzer(self);
        end
end
end
