
function [runtime, tranSysState, tranSync] = runAnalyzer(self, options)
% [RUNTIME, TRANSYSSTATE] = RUN()

T0=tic;
if nargin<2 %~exist('options','var')
    options = self.getOptions;
end
% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
options.timeout_tic = T0;

% options.events (DES event budget) overrides options.samples when set;
% samples remains accepted as a deprecated alias for the event budget.
if isfield(options,'events') && ~isempty(options.events) && isfinite(options.events) && options.events > 0
    options.samples = round(options.events);
end

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'SSA');
if pyHandled
    runtime = pyRuntime;
    return
end

self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);


% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

%options.lang = 'java';

sn = getStruct(self);
% see _kb/06-solver-catalog.md for rationale (SSA FCR)

% Native fork-join support: simulate the tag-augmented copy and fold the
% auxiliary sibling classes back into the original classes at the end
isFJ = any(sn.nodetype == NodeType.Fork) || any(sn.nodetype == NodeType.Join);
if isFJ
    if strcmp(options.lang,'java')
        line_warning(mfilename,'Fork-join models are not supported by the JLINE SSA backend, switching to lang=matlab.\n');
        options.lang = 'matlab';
    end
    if strcmp(options.method,'parallel')
        line_warning(mfilename,'The parallel method does not support fork-join models, switching to the serial method.\n');
    end
    options.method = 'serial';
    sn_orig = sn;
    Korig = sn.nclasses;
    [~, fjsn, fjclassmap] = ModelAdapter.fjtag(self.model);
    sn = fjsn;
    line_debug(options, 'SSA: fork-join tag augmentation, %d classes (%d auxiliary), %d fork firings', sn.nclasses, sn.nclasses-Korig, length(sn.fjsync));
end

switch options.lang
    case 'java'
        line_debug(options, 'SSA: using lang=java, delegating to JLINE');
        % see _kb/06-solver-catalog.md for rationale (SSA lang=java transient trajectory)
        if nargout > 1
            options.method = 'serial';
        end
        switch options.method
            case {'default','serial','parallel'}
                switch options.method
                    case 'default'
                        options.verbose = VerboseLevel.SILENT;
                        actualmethod = 'parallel';
                end
                jmodel = LINE2JLINE(self.model);
                M = jmodel.getNumberOfStations;
                R = jmodel.getNumberOfClasses;
                tic;
                jsolver = JLINE.SolverSSA(jmodel, options);
                [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);

                % see _kb/06-solver-catalog.md for rationale (SSA lang=java transient trajectory)
                if nargout > 1
                    jres = jsolver.result;
                    if isempty(jres.tranSysState)
                        line_error(mfilename,'SSA transient trajectory unavailable from JLINE (serial analyzer required).');
                    end
                    tranSysState = cell(1, sn.nstateful + 1);
                    tranSysState{1} = JLINE.from_jline_matrix(jres.tranSysState.get(java.lang.Integer(0)));
                    for isf = 1:sn.nstateful
                        tranSysState{1+isf} = JLINE.from_jline_matrix(jres.tranSysState.get(java.lang.Integer(isf)));
                    end
                    tranSync = JLINE.from_jline_matrix(jres.tranSync);
                    tranSync = tranSync(:)';
                    % see _kb/06-solver-catalog.md for rationale (SSA lang=java transient trajectory)
                    runtime = jsolver.result.runtime;
                    return;
                end
                CN = JLINE.from_jline_matrix(jsolver.getAvgSysRespT());
                XN = JLINE.from_jline_matrix(jsolver.getAvgSysTput());
                runtime = jsolver.result.runtime;
                QN = reshape(QN',R,M)';
                UN = reshape(UN',R,M)';
                RN = reshape(RN',R,M)';
                TN = reshape(TN',R,M)';
                WN = reshape(WN',R,M)';
                AN = reshape(AN',R,M)';
                % Extract cache hit/miss probabilities from Java model
                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        hitRatioVec = JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getHitRatio());
                        missRatioVec = JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getMissRatio());
                        hitClass = self.model.nodes{ind}.getHitClass;
                        nk = length(hitClass);
                        hitprob = zeros(1, nk);
                        missprob = zeros(1, nk);
                        for k = 1:nk
                            if hitClass(k) > 0 && k <= length(hitRatioVec)
                                hitprob(k) = hitRatioVec(k);
                                missprob(k) = missRatioVec(k);
                            end
                        end
                        self.model.nodes{ind}.setResultHitProb(hitprob);
                        self.model.nodes{ind}.setResultMissProb(missprob);
                    end
                end
                if any(sn.nodetype == NodeType.Cache)
                    self.model.refreshStruct(true);
                end
                self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime);

                % Extract confidence intervals from Java solver results
                if confintEnabled
                    result = jsolver.result;
                    % Helper to check for non-null Java objects
                    isValidMatrix = @(x) ~isempty(x) && isa(x, 'jline.util.matrix.Matrix');
                    if isValidMatrix(result.QNCI)
                        QNCI = JLINE.from_jline_matrix(result.QNCI);
                        QNCI = reshape(QNCI', R, M)';
                    else
                        QNCI = [];
                    end
                    if isValidMatrix(result.UNCI)
                        UNCI = JLINE.from_jline_matrix(result.UNCI);
                        UNCI = reshape(UNCI', R, M)';
                    else
                        UNCI = [];
                    end
                    if isValidMatrix(result.RNCI)
                        RNCI = JLINE.from_jline_matrix(result.RNCI);
                        RNCI = reshape(RNCI', R, M)';
                    else
                        RNCI = [];
                    end
                    if isValidMatrix(result.TNCI)
                        TNCI = JLINE.from_jline_matrix(result.TNCI);
                        TNCI = reshape(TNCI', R, M)';
                    else
                        TNCI = [];
                    end
                    if isValidMatrix(result.ANCI)
                        ANCI = JLINE.from_jline_matrix(result.ANCI);
                        ANCI = reshape(ANCI', R, M)';
                    else
                        ANCI = [];
                    end
                    if isValidMatrix(result.WNCI)
                        WNCI = JLINE.from_jline_matrix(result.WNCI);
                        WNCI = reshape(WNCI', R, M)';
                    else
                        WNCI = [];
                    end
                    % Store CI results
                    self.setAvgResultsCI(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, [], []);
                end
            otherwise
                line_error(mfilename, ['the ',options.method',' method is not available.']);
        end
    case 'matlab'
        line_debug(options, 'SSA: using lang=matlab');
        [QN,UN,RN,TN,CN,XN,~,actualmethod,tranSysState, tranSync, sn, QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = solver_ssa_analyzer(sn, options);

        for isf=1:sn.nstateful
            ind = sn.statefulToNode(isf);
            switch sn.nodetype(sn.statefulToNode(isf))
                case NodeType.Cache
                    self.model.nodes{sn.statefulToNode(isf)}.setResultHitProb(sn.nodeparam{ind}.actualhitprob);
                    self.model.nodes{sn.statefulToNode(isf)}.setResultMissProb(sn.nodeparam{ind}.actualmissprob);
                    if isfield(sn.nodeparam{ind}, 'actualresidt')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultResidT(sn.nodeparam{ind}.actualresidt);
                    end
                    self.model.refreshChains();
            end
        end
        line_debug(options, 'SSA analysis complete: extracting results (nstations=%d, nclasses=%d)', sn.nstations, sn.nclasses);
        runtime = toc(T0);
        T = getAvgTputHandles(self);
        if isFJ
            [QN,UN,RN,TN,CN,XN] = sn_fj_foldback(QN,UN,RN,TN,CN,XN,fjclassmap,Korig);
            self.result.fjclassmap = fjclassmap;
            [TN,~,RN] = sn_pn_avg_rates(sn_orig, QN, TN, [], RN);
            AN = sn_get_arvr_from_tput(sn_orig, TN, T);
            % Join stations report the per-sibling waiting time (JMT
            % convention): QLen over the sibling arrival rate
            for ist=1:sn_orig.nstations
                if sn_orig.nodetype(sn_orig.stationToNode(ist)) == NodeType.Join
                    for r=1:Korig
                        if AN(ist,r) > 0
                            RN(ist,r) = QN(ist,r)/AN(ist,r);
                        end
                    end
                end
            end
        else
            [TN,~,RN] = sn_pn_avg_rates(sn, QN, TN, [], RN);
            AN = sn_get_arvr_from_tput(sn, TN, T);
        end
        if strcmp(options.method,'default') && exist('actualmethod','var')
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod]);
        else
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method);
        end
        self.result.space = sn.space;

        % Store CI data if computed
        if confintEnabled && ~isempty(QNCI)
            self.setAvgResultsCI(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, [], []);
        end
        if lineTimeoutExceeded(options)
            self.result.Avg.timedOut = true;
            line_warning(mfilename,'Solver stopped after the wall-clock time budget (options.timeout=%gs) was exceeded; returning the interim solution.\n', options.timeout);
        end
end
end