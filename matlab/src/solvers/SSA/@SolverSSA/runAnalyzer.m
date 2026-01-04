
function [runtime, tranSysState, tranSync] = runAnalyzer(self, options)
% [RUNTIME, TRANSYSSTATE] = RUN()

T0=tic;
if nargin<2 %~exist('options','var')
    options = self.getOptions;
end
self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

% Show library attribution if verbose and not yet shown
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    sn_temp = self.getStruct();
    libs = SolverSSA.getLibrariesUsed(sn_temp, options);
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

%options.lang = 'java';

sn = getStruct(self);

switch options.lang
    case 'java'
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
                CN = JLINE.from_jline_matrix(jsolver.getAvgSysRespT());
                XN = JLINE.from_jline_matrix(jsolver.getAvgSysTput());
                runtime = jsolver.result.runtime;
                QN = reshape(QN',R,M)';
                UN = reshape(UN',R,M)';
                RN = reshape(RN',R,M)';
                TN = reshape(TN',R,M)';
                WN = reshape(WN',R,M)';
                AN = reshape(AN',R,M)';
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
        [QN,UN,RN,TN,CN,XN,~,actualmethod,tranSysState, tranSync, sn, QNCI, UNCI, RNCI, TNCI, ANCI, WNCI] = solver_ssa_analyzer(sn, options);

        for isf=1:sn.nstateful
            ind = sn.statefulToNode(isf);
            switch sn.nodetype(sn.statefulToNode(isf))
                case NodeType.Cache
                    self.model.nodes{sn.statefulToNode(isf)}.setResultHitProb(sn.nodeparam{ind}.actualhitprob);
                    self.model.nodes{sn.statefulToNode(isf)}.setResultMissProb(sn.nodeparam{ind}.actualmissprob);
                    self.model.refreshChains();
            end
        end
        runtime = toc(T0);
        T = getAvgTputHandles(self);
        AN = sn_get_arvr_from_tput(sn, TN, T);
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
end
end