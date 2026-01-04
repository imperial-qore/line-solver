function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;

if nargin<2
    options = self.getOptions;
end

if ~isinf(options.timespan(1)) && (options.timespan(1) == options.timespan(2))
    line_warning(mfilename,'%s: timespan is a single point, spacing by options.tol (%e).\n',mfilename, options.tol);
    options.timespan(2) = options.timespan(1) + options.tol;
end


self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

% Show library attribution if verbose and not yet shown
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    libs = SolverCTMC.getLibrariesUsed([], options);
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

if self.enableChecks && ~self.supports(self.model)
    line_error(mfilename,'This model contains features not supported by the solver.\n');
end

% Inform user about reducible routing handling
if self.enableChecks
    [isErg, ergInfo] = self.model.isRoutingErgodic();
    if ~isErg && ~isempty(ergInfo.absorbingStations)
        absNames = strjoin(ergInfo.absorbingStations, ', ');
        line_printf([...
            'Note: Model has reducible routing with absorbing stations: %s\n' ...
            '  Results represent limiting/absorption probabilities.\n' ...
            '  Use model.getReducibilityInfo() for detailed analysis.\n'], ...
            absNames);
    end
end

sn = getStruct(self);

% Convert non-Markovian distributions to PH
sn = sn_nonmarkov_toph(sn, options);

M = sn.nstations;
K = sn.nclasses;
NK = sn.njobs;
sizeEstimator = 0;
for k=1:K
    sizeEstimator = sizeEstimator + gammaln(1+NK(k)+M-1) - gammaln(1+M-1) - gammaln(1+NK(k)); % worst-case estimate of the state space
end

if any(isinf(sn.njobs))
    if isinf(options.cutoff)
        line_warning(mfilename,sprintf('The model has open chains, it is recommended to specify a finite cutoff value, e.g., SolverCTMC(model,''cutoff'',1).\n'));
        self.options.cutoff= ceil(6000^(1/(M*K)));
        options.cutoff= ceil(6000^(1/(M*K)));
        line_warning(mfilename,sprintf('Setting cutoff=%d.\n',self.options.cutoff));
    end
    % Mandatory truncation warning for open/mixed models
    line_printf('CTMC solver using state space cutoff = %d for open/mixed model.\n', options.cutoff);
    line_warning(mfilename,'State space truncation may cause inaccurate results. Consider varying cutoff to assess sensitivity.\n');
end

if sizeEstimator > 6
    if ~isfield(options,'force') || options.force == false
        %        line_error(mfilename,'CTMC size may be too large to solve. Stopping SolverCTMC. Set options.force=true to bypass this control.');
        line_error(mfilename,'CTMC size may be too large to solve. Stopping SolverCTMC. Set options.force=true or use SolverCTMC(...,''force'',true) to bypass this control.\n');
        return
    end
end

% we compute all metrics anyway because CTMC has essentially
% the same cost
if isinf(options.timespan(1))
    s0 = sn.state;
    s0prior = sn.stateprior;
    for ind=1:sn.nnodes
        if sn.isstateful(ind)
            isf = sn.nodeToStateful(ind);
            sn.state{isf} = s0{isf}(maxpos(s0prior{1}),:); % pick one particular initial state
        end
    end
    [QN,UN,RN,TN,CN,XN,Q,SS,SSq,Dfilt,~,~,sn] = solver_ctmc_analyzer(sn, options);
    % update initial state if this has been corrected by the state space
    % generator
    for isf=1:sn.nstateful
        ind = sn.statefulToNode(isf);
        self.model.nodes{ind}.setState(sn.state{isf});
        switch class(self.model.nodes{sn.statefulToNode(isf)})
            case 'Cache'
                self.model.nodes{sn.statefulToNode(isf)}.setResultHitProb(sn.nodeparam{ind}.actualhitprob);
                self.model.nodes{sn.statefulToNode(isf)}.setResultMissProb(sn.nodeparam{ind}.actualmissprob);
                self.model.refreshChains();
        end
    end
    %sn.space = SS;
    self.result.infGen = Q;
    self.result.space = SS;
    self.result.spaceAggr = SSq;
    self.result.nodeSpace = sn.space;
    self.result.eventFilt = Dfilt;
    runtime = toc(T0);
    sn.space = {};
    T = getAvgTputHandles(self);
    AN=sn_get_arvr_from_tput(sn, TN, T);
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method);
else
    lastSol= [];
    s0 = sn.space;
    s0prior = sn.stateprior;

    s0_sz = cellfun(@(x) size(x,1), s0)';
    s0_id = pprod(s0_sz-1);
    cur_state = sn.state;
    while s0_id>=0 % for all possible initial states
        s0prior_val = 1;
        for ind=1:sn.nnodes
            if sn.isstateful(ind)
                isf = sn.nodeToStateful(ind);
                s0prior_val = s0prior_val * s0prior{isf}(1+s0_id(isf)); % update prior
                sn.state{isf} = s0{isf}(1+s0_id(isf),:); % assign initial state to network
            end
        end
        if s0prior_val > 0
            [t,pit,QNt,UNt,~,TNt,~,~,Q,SS,SSq,Dfilt,runtime_t] = solver_ctmc_transient_analyzer(sn, options);
            self.result.space = SS;
            self.result.spaceAggr = SSq;
            self.result.infGen = Q;
            self.result.eventFilt = Dfilt;
            %sn.space = SS;
            setTranProb(self,t,pit,SS,runtime_t);
            if isempty(self.result) || ~isfield(self.result,'Tran') || ~isfield(self.result.Tran,'Avg') || ~isfield(self.result.Tran.Avg,'Q')
                self.result.Tran.Avg.Q = cell(M,K);
                self.result.Tran.Avg.U = cell(M,K);
                self.result.Tran.Avg.T = cell(M,K);
                for ist=1:M
                    for r=1:K
                        self.result.Tran.Avg.Q{ist,r} = [QNt{ist,r} * s0prior_val,t];
                        self.result.Tran.Avg.U{ist,r} = [UNt{ist,r} * s0prior_val,t];
                        self.result.Tran.Avg.T{ist,r} = [TNt{ist,r} * s0prior_val,t];
                    end
                end
            else
                for ist=1:M
                    for r=1:K
                        tunion = union(self.result.Tran.Avg.Q{ist,r}(:,2), t);
                        dataOld = interp1(self.result.Tran.Avg.Q{ist,r}(:,2),self.result.Tran.Avg.Q{ist,r}(:,1),tunion);
                        dataNew = interp1(t,QNt{ist,r},tunion);
                        self.result.Tran.Avg.Q{ist,r} = [dataOld+s0prior_val*dataNew,tunion];
                        dataOld = interp1(self.result.Tran.Avg.U{ist,r}(:,2),self.result.Tran.Avg.U{ist,r}(:,1),tunion);
                        dataNew = interp1(t,UNt{ist,r},tunion);
                        self.result.Tran.Avg.U{ist,r} = [dataOld+s0prior_val*dataNew,tunion];

                        dataOld = interp1(self.result.Tran.Avg.T{ist,r}(:,2),self.result.Tran.Avg.T{ist,r}(:,1),tunion);
                        dataNew = interp1(t,TNt{ist,r},tunion);
                        self.result.Tran.Avg.T{ist,r} = [dataOld+s0prior_val*dataNew,tunion];
                    end
                end
            end
        end
        s0_id=pprod(s0_id,s0_sz-1); % update initial state
    end
    % Now we restore the original state
    for ind=1:sn.nnodes
        if sn.isstateful(ind)
            isf = sn.nodeToStateful(ind);
            self.model.nodes{ind}.setState(cur_state{isf});
        end
    end

    runtime = toc(T0);
    sn.space = {};
    self.result.('solver') = getName(self);
    self.result.runtime = runtime;
    self.result.solverSpecific = lastSol;
end
end