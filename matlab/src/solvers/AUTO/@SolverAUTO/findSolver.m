function T = findSolver(self, metric, showAll)
% T = FINDSOLVER(METRIC, SHOWALL)
%
% Which solvers and solver methods can analyze this model, and for the ones
% that cannot, why not.
%
% One row per (family, method) pair AUTO can be asked for, with columns
%
%   Solver    the method family, 'mva', 'ctmc', 'ldes', ...
%   Method    the method name to pass as a method name, 'mva.exact'
%   Runnable  true when the model passes that method's own support gate
%   Class     'exact', 'approx', 'bound' or 'simulation'
%   Metrics   the measure groups the family answers, see METRICGROUPS
%   Reason    why a refused pair was refused, '' when RUNNABLE
%
% METRIC narrows the report to the pairs that answer one measure, named either
% by its group ('cdf') or by the accessor that returns it ('getCdfRespT'). ''
% or 'any' keeps every pair. SHOWALL keeps the refused pairs too; by default
% only the runnable ones are listed, since a caller asking what it can run has
% no use for the two hundred rows that say it cannot.
%
% THE GATE IS NOT A SECOND ONE. It is the gate CHOOSESOLVERRANKED applies
% before delegating, asked of every candidate instead of of the first feasible
% one, which is exactly what LISTVALIDMETHODS already does -- that method is
% now the Method column of the runnable rows, so the two cannot disagree.
% What is new is that the REASON the gate produced is kept rather than
% discarded, and that the answer carries the two facts a caller needs in order
% to choose among the survivors: whether the method is exact on this model,
% and which measures it can report.
%
% AND IT IS THE RUN'S GATE, ASKED OF THE RESOLVED NAME. A name is gated the
% way NetworkSolver.runAnalyzerChecks gates it, through
% SUPPORTSRESOLVEDMETHOD: 'default' is first resolved to the method it runs
% as on this model (rqna on a MAP-fed open queue in SolverMVA, mem on a
% non-Markovian open model in SolverNC) and the gate is asked of that. Asking
% the literal name against the base envelope made the report and the run
% disagree in both directions: mva.default was refused on a MAP-fed queue that
% ran through rqna, and offered on a MAP-fed fork-join that the run, resolving
% to rqna, refused. When the two names differ, the Reason of a refused row
% says which method the name ran as.
%
% WHY A REASON HAS TO BE RECONSTRUCTED HERE for some rows. The base
% SUPPORTSMODELMETHOD returns a reason only when the solver diverges per
% method; when it does not, it falls back to SUPPORTS(MODEL), which answers
% with a bare logical. That is enough for a gate, which only has to stop the
% run, and not enough for a report, whose whole content is the explanation. So
% a refused row is re-asked against the solver's own feature set, which is
% where the offending feature names are, and THE FEATURE-SET REASON WINS
% WHENEVER THERE IS ONE: a gate may answer with a sentence that names no
% feature, and that does not tell a user what to change. A refusal the feature
% set does NOT explain keeps the gate's own words, because then the gate knows
% something the feature set cannot express (product form, a binding buffer,
% NC 'mem' applicability).
%
% A REPORT MUST NOT PRINT. Asking a solver whether it supports the model runs
% SOLVERFEATURESET.SUPPORTS, which warns naming the missing feature: a side
% effect that is right on the solve path, where nobody asked to be told, and
% wrong here, where every refused row would raise one and the answer IS the
% table. The guard below silences the walk and restores the caller's level on
% every exit path, error included.
%
% See also SolverAUTO.listValidMethods, Network.findSolver, SolverFeatureSet
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    metric = '';
end
if nargin < 3
    showAll = false;
end

% Normalized to char FIRST: a string scalar "" is not empty by ISEMPTY, so a
% findSolver("") that means "every measure" would otherwise fall through to the
% unknown-measure error below.
if isstring(metric)
    metric = char(metric);
end
group = SolverAUTO.metricGroupOf(metric);
if isempty(group) && ~isempty(metric) && ~any(strcmpi(metric,{'any','all'}))
    line_error(mfilename, sprintf(['''%s'' names no measure. Pass a group (%s) ' ...
        'or the accessor that returns it, e.g. ''getCdfRespT''.'], ...
        metric, strjoin(SolverAUTO.metricGroups(), ', ')));
end

% The three model properties an exactness claim can rest on, evaluated once: a
% method whose exactness needs one of them reports 'approx' without it, see
% SolverAUTO.methodClass.
isProductForm = false;
isQbdShape = false;
hasCache = false;
if isa(self.model,'Network')
    isProductForm = self.model.hasProductFormSolution();
    sn = self.model.getStruct(false);
    nsources = sum(sn.nodetype == NodeType.Source);
    isQbdShape = all(isinf(sn.njobs)) && (sn.nstations - nsources) == 1;
    hasCache = any(sn.nodetype == NodeType.Cache);
end

verboseGuard = GlobalConstants.pushVerbose(VerboseLevel.SILENT); %#ok<NASGU>

% Column-shaped from the start, so that a model no family serves still returns
% a table with the six variables rather than one TABLE could not build: a 0x0
% cell is not a zero-row column, and an empty report has to be a table a caller
% can index and concatenate like any other.
Solver = cell(0,1);
Method = cell(0,1);
Runnable = false(0,1);
Class = cell(0,1);
Metrics = cell(0,1);
Reason = cell(0,1);

probeOptions = self.options;
probeOptions.verbose = 0;
probeOptions.method = 'default';

families = SolverAUTO.familyNames();
for f = 1:length(families)
    family = families{f};
    if ~SolverAUTO.familyAcceptsModelClass(family, self.model)
        continue
    end
    groups = SolverAUTO.familyMetrics(family);
    if ~isempty(group) && ~any(strcmp(group, groups))
        continue
    end
    try
        probe = self.buildFamilySolver(family, self.model, probeOptions);
        declared = probe.listValidMethods();
    catch
        % A family that cannot even be instantiated here (SolverLQNS without
        % the lqns binary, SolverLN on a flat Network) contributes nothing:
        % there is no solver to report on and no gate to ask.
        continue
    end
    metricList = strjoin(groups, ',');
    % Once per family, not once per refused row: the flat feature set does not
    % vary with the method, and a family that refuses every one of its forty
    % methods would otherwise recompute the same comparison forty times.
    famFeatureReason = featureReason(probe, self.model);
    for m = 1:length(declared)
        name = declared{m};
        if length(name) > length(family) && strncmp(name, [family,'.'], length(family)+1)
            % A spelling already qualified with its own family. SolverFLD
            % declares both 'dae' and 'fluid.dae' so that its own gate takes
            % either, and prefixing the family again yields 'fluid.fluid.dae':
            % a method name that does resolve, but that names the same method twice
            % and would double every fluid row of this report.
            continue
        end
        if SolverAUTO.isMethodAlias(family, name, declared)
            % The same duplication under a different prefix. SolverMVA
            % advertises every AMVA name twice, plain and 'amva.'-prefixed,
            % and its dispatch strips the prefix, so the two spellings are one
            % algorithm; that alone was 20 of the 49 mva rows of a report. The
            % plain spelling is the one kept.
            continue
        end
        [ok, reason, resolved] = gateVerdict(probe, self.model, name);
        if ok
            reason = '';
        elseif isGenericReason(reason)
            % A SPECIFIC REASON WINS; the feature-set names only replace a
            % generic one or fill an empty one. Getting this the other way
            % round was worse than saying nothing: UQ's feature set answers a
            % DIFFERENT question -- it declares the Prior that UQ itself
            % consumes, every other feature being the INNER solver's to accept
            % (see UQ.supports) -- so comparing it against the model listed
            % every feature the model uses and buried the real reason, that
            % the model carries no uncertain parameter at all.
            if ~isempty(famFeatureReason)
                reason = famFeatureReason;
            elseif isempty(reason)
                reason = 'This solver refuses the model through its own structural check.';
            end
        end
        if ~strcmp(resolved, name)
            % The verdict is about the method the name RESOLVED to, which a
            % reader of the row cannot otherwise know: 'default' refused for
            % "(feature: Join)" reads as nonsense until it says it ran as rqna,
            % and a RUNNABLE 'mfq' on a tandem is honest only if the row says
            % it runs as matrix (the documented fallback, which the run keeps).
            reason = strtrim(sprintf('''%s'' runs as ''%s'' on this model. %s', name, resolved, reason));
        end
        if ~ok && ~showAll
            continue
        end
        Solver{end+1,1} = family; %#ok<AGROW>
        Method{end+1,1} = [family,'.',name]; %#ok<AGROW>
        Runnable(end+1,1) = ok; %#ok<AGROW>
        Class{end+1,1} = SolverAUTO.methodClass(family, name, ...
            stochasticVerdict(probe, name), isProductForm, isQbdShape, hasCache); %#ok<AGROW>
        Metrics{end+1,1} = metricList; %#ok<AGROW>
        Reason{end+1,1} = reason; %#ok<AGROW>
    end
end

T = table(Solver, Method, Runnable, Class, Metrics, Reason, ...
    'VariableNames', {'Solver','Method','Runnable','Class','Metrics','Reason'});
end

function bool = isGenericReason(reason)
% BOOL = ISGENERICREASON(REASON)
%
% Does REASON say only that SOME feature is unsupported, without naming one?
% That is the sentence SolverNC and the shared base gate return, and it is the
% one worth replacing with the offending feature names. Empty counts: a gate
% that answered with a bare logical said nothing at all.
bool = isempty(reason) || (contains(lower(reason), 'features are not supported') && ...
    ~contains(reason, '(feature:'));
end

function [ok, reason, resolved] = gateVerdict(probe, model, method)
% [OK, REASON, RESOLVED] = GATEVERDICT(PROBE, MODEL, METHOD)
%
% The method-level gate, asked without letting it raise. The base
% SUPPORTSMODELMETHOD already falls back to the flat SUPPORTS(MODEL) for a
% solver that does not diverge per method, so this one call is both gates.
%
% RESOLVED is the method METHOD runs as on this model. A NetworkSolver is
% asked through SUPPORTSRESOLVEDMETHOD, the very call RUNANALYZERCHECKS makes,
% so that 'default' is gated as the method it resolves to and the report
% cannot call Runnable a pair the run refuses, or refuse one it serves. The
% probe carries the options the run would carry (AUTO's own, method reset),
% which is what a resolution may read.
%
% MODEL is passed in rather than read off the probe: a LayeredNetwork solver
% keeps its model under another name, and SUPPORTS is a question about the
% model this report is about in any case.
ok = true;
reason = '';
resolved = method;
if ~ismethod(probe, 'supportsModelMethod')
    try
        % A LayeredNetwork solver answers [bool, reason]; take the first
        % output, as chooseSolverRanked does.
        ok = probe.supports(model);
    catch
        ok = true; % no claim, no gate
    end
    return
end
try
    if ismethod(probe, 'supportsResolvedMethod')
        opts = probe.getOptions();
        opts.method = method;
        [ok, reason, resolved] = probe.supportsResolvedMethod(opts);
    else
        % An ensemble solver (LN, UQ) has no resolution step: its name is the
        % method that runs.
        [ok, reason] = probe.supportsModelMethod(method);
    end
catch ME
    % A gate that raises has said something, and it is the only thing it can
    % say about this pair; reporting it beats swallowing it and calling the
    % pair runnable.
    ok = false;
    reason = ME.message;
end
end

function bool = stochasticVerdict(probe, method)
% BOOL = STOCHASTICVERDICT(PROBE, METHOD)
% ISSTOCHASTICMETHOD asked without letting it raise; a solver that cannot
% classify a name is taken at its class default, deterministic.
bool = false;
try
    bool = probe.isStochasticMethod(method);
catch
end
end

function reason = featureReason(probe, model)
% REASON = FEATUREREASON(PROBE, MODEL)
%
% The offending feature names, or '' when the solver's feature set accepts the
% model.
%
% The empty answer is meaningful and not a failure: it says the refusal came
% from somewhere the feature set cannot see, so the caller keeps whatever the
% gate itself said. The two outputs of SUPPORTS are both taken, which is also
% what suppresses its side-effect warning.
reason = '';
mc = metaclass(probe);
if ~any(strcmp({mc.MethodList.Name},'getFeatureSet'))
    return
end
featUsed = [];
try
    featUsed = model.getUsedLangFeatures();
catch
end
if ~isa(featUsed, 'SolverFeatureSet')
    % A LayeredNetwork answers with one feature set PER LAYER, which is not a
    % set this comparison can take; the gate's own words stand instead.
    return
end
try
    featSupported = feval([class(probe),'.getFeatureSet']);
catch
    return
end
% A SOLVER THAT DECLARES NOTHING HAS NO ENVELOPE, and a missing-feature list
% against an empty set is not an explanation: it is every feature the model
% uses. UQ is the case -- it computes nothing itself, it expands a Prior and
% runs another solver at each design point -- and it refused an M/M/1 with
% "Some features are not supported (feature: Sink, Source, Exp,
% RoutingStrategy_PROB, SchedStrategy_FCFS, OpenClass)", which tells a user
% nothing they can act on. Such a refusal is structural, so the caller keeps
% the gate's own words or the structural sentence.
declared = struct2cell(featSupported.list);
if ~any([declared{:}])
    return
end
try
    [~, reason] = SolverFeatureSet.supports(featSupported, featUsed);
catch
    reason = '';
end
end
