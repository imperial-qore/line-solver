function [runtime, analyzer] = runAnalyzer(self, options)
% RUNANALYZER Execute the bound-analysis solver.
%
% Dispatches noniterative bound families to the shared bound handler
% solver_ba_analyzer, and hierarchical/iterative families to their
% SolverBA-native analyzers. Each method returns a single (upper or lower)
% bound; use getBounds() to obtain the {lower,upper} bracket for a family.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = self.getOptions;
end

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'BA');
if pyHandled
    runtime = pyRuntime;
    analyzer = [];
    return
end

% Solver console: SolverBA does not pass through runAnalyzerChecks, so it
% opens its own run here. The guard must live until this function returns.
consoleGuard = LineConsole.beginRun(self, options); %#ok<NASGU>

sn = self.getStruct;

% Every bound family here is parameterized by demands and a closed population,
% except the three OPEN-network families 'bpt', 'bgt' and 'snc', which are
% derived for an open network and are refused on a closed one instead, and the
% 'spnlp' family, which is parameterized by a MARKING: whether that marking is
% bounded is a question about the P-invariants of the net and not about
% sn.nclosedjobs, so spn_lpbnd decides it rather than this gate.
if sn.nclosedjobs <= 0 && ~startsWith(options.method,'bpt') && ...
        ~startsWith(options.method,'bgt') && ~startsWith(options.method,'snc') && ...
        ~startsWith(options.method,'spnlp')
    line_error(mfilename, 'SolverBA supports closed queueing networks only.');
end

% 'default' -> geometric upper bound; 'qr' -> QRF quadratic reduction; 'lr' ->
% LP linear reduction (bare 'lr' means 'lr.upper'), distinct from
% 'qrf.mmi.linear'. On a model with a binding finite buffer, 'default',
% 'auto' and 'auto.upper' are routed to 'qrf.bas' where that bound applies
% (the same routing SolverMVA does with sn_is_bas_model -> 'sqd'). Resolved
% through the shared BA_RESOLVE_MODEL_METHOD so that resolveMethod (hence
% findSolver), listValidMethods, getMethodFeatureSet and the structural
% predicate all gate on the same name this dispatches on.
% see _kb/06-solver-catalog.md for rationale
method = ba_resolve_model_method(sn, options.method);
if ~strcmp(method, ba_resolve_method(options.method))
    line_debug(options, sprintf(['BA: %s resolves to ''%s'' on this model, which has ' ...
        'a binding finite buffer'], options.method, method));
end

% THE GATE, THE SAME ONE MODEL.HELP REPORTS. SolverBA does not pass through
% runAnalyzerChecks (that gate refuses a name outside the MODEL-NARROWED
% listValidMethods, whereas a direct request must reach the analyzer's own
% reason), so supportsModelMethod is asked here directly, with its two halves:
% the feature envelope of the resolved method (getMethodFeatureSet, through
% the base gate) and the structural predicate BA_METHOD_REFUSAL. Until this
% was here the envelope was report-only -- a MAP-serviced or fork-join model
% refused by model.help was still bounded on its means when asked by name --
% and the blocking, product-form and routing rules lived in the analyzers.
% Placed after the aliases so 'default' is judged as what it runs as, and
% before the lang switch so lang='java' is refused on the same terms.
% supportsModelMethod asks the structural predicate FIRST and the feature
% envelope second, so a model outside both gates is refused with the
% analyzer's own sentence ("supports fully open networks only") rather than
% with a feature list. Its FORRUN flag skips BA_METHOD_DEGENERATE: a vacuous
% bound is a valid one, withheld from the report but published when named.
if self.enableChecks
    [ok, reason] = self.supportsModelMethod(method, true);
    if ~ok
        if strcmp(method, options.method)
            line_error(mfilename, sprintf('This model contains features not supported by the solver. %s', reason));
        else
            line_error(mfilename, sprintf('This model contains features not supported by the solver''s ''%s'' method. %s', method, reason));
        end
    end
end

% Show library attribution (QRF uses the Optimization Toolbox) once.
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    libs = SolverBA.getLibrariesUsed([], setfield(options,'method',method)); %#ok<SFLD>
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

% lang='java' delegates the whole bound solve to the JAR, which is how a MATLAB
% user reaches jline.api.pfqn.mva.Pfqn_ssd: Solver_mva_bound_analyzer routes
% every method starting with 'ssd' to it. Guard placed here, after the method
% aliases are resolved and before any MATLAB analyzer runs, matching
% @SolverNC/runAnalyzer.m so there is one convention rather than two.
switch options.lang
    case 'java'
        line_debug(options, 'BA: using lang=java, delegating to JLINE');
        jmodel = LINE2JLINE(self.model);
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jopts = options; jopts.method = method;
        jsolver = JLINE.SolverBA(jmodel, jopts);
        % getAvgTable(true) is the UNFILTERED grid, and the reshape below needs it:
        % the no-argument getter DROPS every (station,class) cell whose six metrics
        % are all zero, so on a model with a disabled pair it returns fewer than M*R
        % entries and reshape(...,R,M) errors out. MATLAB applies its own filter when
        % the table is PRINTED, so the bridge must carry the whole grid, zeros included.
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable(true));
        runtime = jsolver.result.runtime;
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        TN = reshape(TN',R,M)';
        WN = reshape(WN',R,M)';
        AN = reshape(AN',R,M)';
        % C and X are per-chain scalars the bound families fill in MATLAB;
        % the JAR table does not carry them, so they stay empty rather than
        % being invented here.
        CN = [];
        XN = [];
        analyzer = [];
        self.model.refreshStruct(true);
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,method,1);
        return
    case 'matlab'
        line_debug(options, 'BA: using lang=matlab');
end

iter = 1;
if startsWith(method,'qrf')
    % QRF (Quadratic Reduction Framework) LP-based bounds for single-class
    % closed networks with PH service. The adapter enforces its own gating.
    analyzer = @(qn) solver_ba_qrf_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,runtime] = solver_ba_qrf_analyzer(sn, bopts);
elseif startsWith(method,'spnlp')
    % Moment-relaxation LP bounds (Liu 1998) for a stochastic Petri net. The
    % adapter enforces its own gating.
    analyzer = @(qn) solver_ba_spnlp_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_ba_spnlp_analyzer(sn, bopts); %#ok<ASGLU>
elseif startsWith(method,'bgt')
    % Piecewise-linear Lyapunov upper bound (Bertsimas-Gamarnik-Tsitsiklis
    % 2001) for multitype OPEN networks. The adapter enforces its own gating.
    analyzer = @(qn) solver_ba_bgt_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_ba_bgt_analyzer(sn, bopts); %#ok<ASGLU>
elseif startsWith(method,'snc')
    % Stochastic network calculus upper bound (Fidler-Rizk 2015) for
    % feed-forward OPEN networks. The adapter enforces its own gating.
    analyzer = @(qn) solver_ba_snc_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_ba_snc_analyzer(sn, bopts); %#ok<ASGLU>
elseif startsWith(method,'bpt')
    % Achievable-region LP relaxation (Bertsimas-Paschalidis-Tsitsiklis 1994)
    % for multiclass OPEN networks. The adapter enforces its own gating.
    analyzer = @(qn) solver_ba_bpt_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,lG,runtime] = solver_ba_bpt_analyzer(sn, bopts); %#ok<ASGLU>
elseif any(strcmp(method, SolverBA.listAllMethods()))
    % listAllMethods, NOT listValidMethods: the latter drops what this model
    % cannot run, and dispatching on it would answer a direct request for a
    % gated method with 'unknown method' instead of the analyzer's own reason.
    analyzer = @(qn) solver_ba_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,lG,runtime,iter] = solver_ba_analyzer(sn, bopts);
else
    line_error(mfilename, ['Unknown bound method ''%s''. Valid methods: %s'], ...
        method, strjoin(self.listValidMethods(), ', '));
end

M = sn.nstations;
R = sn.nclasses;
% Arrival rates from throughputs via sn_get_arvr_from_tput (not a zero matrix).
% see _kb/06-solver-catalog.md for rationale
AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
WN = zeros(M,R);
self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,method,iter);
end
