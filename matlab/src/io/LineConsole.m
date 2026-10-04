classdef LineConsole
    % LINECONSOLE Running progress log of a LINE solver run.
    %
    % LineConsole narrates what a solver is doing while it does it: reading
    % the model, compiling the network structure, computing chains, visits and
    % demands, resolving the method, iterating, and closing with the figures
    % of merit. Each line carries the elapsed time since the run started:
    %
    %   [   0.014s] compiling the network structure (sn)
    %   [   0.031s]   computing chains and visit ratios
    %   [   0.052s] closed queueing network: 3 stations, 1 class, 1 chain
    %   [   0.061s]   AMVA sweep 10: residual 1.11e-01, X = 1.2225
    %
    % It prints no tables: the result table stays the caller's own
    % getAvgTable. THE CONSOLE IS VerboseLevel.DEBUG: it narrates exactly when
    % the run is at DEBUG and is silent at every lower level, and it never
    % alters a numerical result. There is no separate console switch -- the
    % console was one, until it became clear that a running progress log IS
    % what a debug verbosity is for, and two switches for one channel only let
    % a session ask for DEBUG and get nothing.
    %
    % Typical use:
    %   line_verbosity(VerboseLevel.DEBUG)
    %   SolverMVA(model).getAvgTable
    % or, for one run alone:
    %   SolverMVA(model,'verbose',VerboseLevel.DEBUG).getAvgTable
    %
    % Nested runs (an inner solver invoked by SolverLN or SolverAUTO) do not
    % narrate: only the outermost run writes, so several layer solves cannot
    % interleave their lines. Use ISACTIVE to suppress a legacy print, OWNSLOG
    % to emit.
    %
    % See also: line_verbosity, VerboseLevel, GlobalConstants
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)

        %% ---------------- state ----------------

        function st = state(newst)
            % ST = STATE(NEWST) reads, and optionally replaces, console state
            persistent CONSOLESTATE
            if isempty(CONSOLESTATE)
                CONSOLESTATE = LineConsole.blankState();
            end
            if nargin > 0
                CONSOLESTATE = newst;
            end
            st = CONSOLESTATE;
        end

        function st = blankState()
            % ST = BLANKSTATE() initial console state
            st.depth = 0;      % nesting level of open solver runs
            st.muted = 0;      % >0 while a silenced run is executing
            st.active = false; % true while an outermost run is narrating
            st.tag = '';       % solver name of the open run
            st.modelName = ''; % name of the model under study
            st.t0 = [];        % tic marker of the open run
            st.tsetup = NaN;   % seconds spent before the analysis began
            st.quiet = 0;      % >0 while detail lines are being suppressed
            st.forceDetail = false; % true while the run compiles its own model
            st.lastLoop = '';  % text of the iteration loop last announced
            st.shown = 0;      % iteration lines written for the current loop
            st.trunc = false;  % true once that loop's lines were cut short
            st.pending = {};   % legacy completion lines held until the run closes
        end

        function reset()
            % RESET() forgets any open run (used after an interrupted solve)
            LineConsole.state(LineConsole.blankState());
        end

        function tf = isActive()
            % TF = ISACTIVE() true while a run is narrating
            st = LineConsole.state();
            tf = st.active;
        end

        function deferPrint(fmt, varargin)
            % DEFERPRINT(FMT, ...) queues a legacy line for after the run closes
            %
            % The console owns the log while a run narrates, so a solver's
            % standard completion message is held here and written once the
            % closing DONE line has gone out, reading as it would with the
            % console off.
            st = LineConsole.state();
            if ~st.active
                line_printf(fmt, varargin{:});
                return
            end
            if st.depth > 1
                return % a nested run does not narrate, and does not report
            end
            st.pending{end+1} = sprintf(fmt, varargin{:});
            LineConsole.state(st);
        end

        function unmute()
            % UNMUTE() ends one muted scope opened by a silenced run
            st = LineConsole.state();
            st.muted = max(0, st.muted - 1);
            LineConsole.state(st);
        end

        function tf = writes()
            % TF = WRITES() true when a progress line should be printed
            %
            % Either the outermost run is narrating, or no run is open and
            % the session is at DEBUG -- the second case is what lets model
            % construction (refreshStruct) narrate before any solver exists.
            st = LineConsole.state();
            if st.muted > 0
                tf = false;
            elseif st.depth > 0
                tf = LineConsole.ownsLog();
            else
                tf = LineConsole.wanted();
            end
        end

        function guard = pushQuiet()
            % GUARD = PUSHQUIET() suppresses detail lines for a bounded scope
            %
            % Used where many auxiliary models are compiled in a row (the
            % SolverLN layer builders): each keeps its one headline STEP and
            % drops its SUBSTEP breakdown. The caller must hold the handle.
            st = LineConsole.state();
            st.quiet = st.quiet + 1;
            LineConsole.state(st);
            guard = onCleanup(@() LineConsole.popQuiet());
        end

        function popQuiet()
            % POPQUIET() ends one PUSHQUIET scope
            st = LineConsole.state();
            st.quiet = max(0, st.quiet - 1);
            LineConsole.state(st);
        end

        function tf = ownsLog()
            % TF = OWNSLOG() true only inside the OUTERMOST open run
            %
            % Analyzer hooks test this rather than ISACTIVE: an inner solver
            % driven by SolverLN or SolverAUTO runs its own loops, and letting
            % it write would interleave several narrations. Use ISACTIVE to
            % suppress a legacy print, OWNSLOG to emit.
            st = LineConsole.state();
            tf = st.active && st.depth <= 1 && st.muted == 0;
        end

        function tf = wanted(options)
            % TF = WANTED(OPTIONS) resolves whether this run should narrate
            %
            % THE CONSOLE IS DEBUG, and this is the whole rule: a run narrates
            % when it is at VerboseLevel.DEBUG and at no lower level. The run's
            % own options.verbose decides when it carries one, otherwise the
            % session level does; there is no third switch that could put the
            % two out of step.
            tf = LineConsole.isDebug(GlobalConstants.Verbose);
            if nargin > 0 && isstruct(options) && isfield(options,'verbose') && ...
                    ~isempty(options.verbose)
                if islogical(options.verbose)
                    % A LOGICAL CANNOT NAME DEBUG. 'verbose',true is
                    % long-standing usage for "not silent" and carries no third
                    % state, so it may VETO the console but never switch it on:
                    % false silences the run, true defers to the session level.
                    % Reading it numerically instead would make true == 1 == STD
                    % pass a `>= DEBUG` test only by accident of the encoding.
                    tf = tf && all(options.verbose(:));
                else
                    tf = LineConsole.isDebug(options.verbose);
                end
            end
        end

        function tf = isDebug(level)
            % TF = ISDEBUG(LEVEL) true when a verbosity LEVEL asks for DEBUG
            %
            % A logical is not a level (see WANTED) and never reaches here from
            % that path; treated as false so a stray one cannot switch the
            % console on by itself.
            if isempty(level) || islogical(level)
                tf = false;
            else
                tf = any(double(level(:)) >= VerboseLevel.DEBUG);
            end
        end

        %% ---------------- run lifecycle ----------------

        function guard = beginRun(self, options)
            % GUARD = BEGINRUN(SELF, OPTIONS) opens a console run
            %
            % Returns an onCleanup handle that closes the run, writing the
            % closing lines, when the caller's analyzer returns or errors. The
            % caller MUST hold the handle for the whole analysis.
            st = LineConsole.state();
            if st.depth == 0 && ~LineConsole.wanted(options)
                % A run that must stay silent MUTES the console for its whole
                % duration, so that nothing it triggers -- a structure compile,
                % a step line of its own -- leaks out at depth 0. SolverLQNS
                % and SolverENV default to a silent verbosity and take this
                % path unless the caller asks for VerboseLevel.DEBUG.
                st.muted = st.muted + 1;
                LineConsole.state(st);
                guard = onCleanup(@() LineConsole.unmute());
                return
            end
            st.depth = st.depth + 1;
            if st.depth > 1 % nested run: the outer analyzer owns the log
                LineConsole.state(st);
                guard = onCleanup(@() LineConsole.closeRun());
                return
            end
            st.active = true;
            st.tag = LineConsole.solverTag(self);
            st.modelName = LineConsole.modelName(self);
            st.t0 = tic;
            st.tsetup = NaN;
            st.shown = 0;
            st.trunc = false;
            LineConsole.state(st);
            LineConsole.detail('', true); % fresh reporting budget per run

            LineConsole.openingLines(self, options);
            st = LineConsole.state();
            st.tsetup = toc(st.t0);
            LineConsole.state(st);
            guard = onCleanup(@() LineConsole.closeRun(self));
        end

        function closeRun(self)
            % CLOSERUN(SELF) closes the innermost open run
            st = LineConsole.state();
            if st.depth <= 0
                return
            end
            st.depth = st.depth - 1;
            LineConsole.state(st);
            if st.depth > 0 || ~st.active
                return
            end
            if nargin > 0
                LineConsole.closingLines(self);
            end
            % the standard completion message follows the closing DONE line
            st = LineConsole.state();
            for p = 1:numel(st.pending)
                line_printf('%s', st.pending{p});
            end
            % the mute of an enclosing silenced run outlives this run's state
            fresh = LineConsole.blankState();
            fresh.muted = st.muted;
            LineConsole.state(fresh);
        end

        %% ---------------- the opening narration ----------------

        function openingLines(self, options)
            % OPENINGLINES(SELF, OPTIONS) narrates the setup of the run
            st = LineConsole.state();
            LineConsole.sessionClock(true); % each run's timeline starts at zero
            if LineConsole.writes() % the opening row is set off from what preceded it
                fprintf(GlobalConstants.StdOut, '\n');
            end
            LineConsole.step('LINE %s: Solver%s starting on model ''%s'' (lang %s)', ...
                LineConsole.versionString(), st.tag, LineConsole.modelName(self), ...
                options.lang);
            LineConsole.readModel(self);
            LineConsole.compileStruct(self);
            LineConsole.recognizeModel(self);
            LineConsole.reportMethod(self, options);
            LineConsole.presolve(self);
        end

        function readModel(self)
            % READMODEL(SELF) what the model object declares, before compiling
            model = LineConsole.modelOf(self);
            if isempty(model) || ~ismethod(model,'getNumberOfNodes')
                return
            end
            LineConsole.step('reading the model: %s, %s', ...
                LineConsole.plural(model.getNumberOfNodes(),'node'), ...
                LineConsole.plural(model.getNumberOfClasses(),'job class','job classes'));
        end

        function compileStruct(self)
            % COMPILESTRUCT(SELF) compiles sn, narrated by refreshStruct
            model = LineConsole.modelOf(self);
            if isempty(model) || ~ismethod(model,'getStruct')
                return
            end
            if LineConsole.hasCompiledStruct(model)
                LineConsole.step('network structure already compiled, reusing it');
                return
            end
            st = LineConsole.state();
            st.forceDetail = true; % this compile is the run's own model
            LineConsole.state(st);
            model.getStruct(); % refreshStruct emits its own substeps
            st = LineConsole.state();
            st.forceDetail = false;
            LineConsole.state(st);
        end

        function tf = hasCompiledStruct(model)
            % TF = HASCOMPILEDSTRUCT(MODEL) true when sn is already built
            tf = isprop(model,'sn') && ~isempty(model.sn);
            if tf && isprop(model,'hasStruct')
                tf = logical(model.hasStruct);
            end
        end

        function recognizeModel(self)
            % RECOGNIZEMODEL(SELF) states what kind of model this is
            sn = LineConsole.structOf(self);
            if isempty(sn)
                return
            end
            if ~LineConsole.has(sn,'nstations') % layered
                if LineConsole.has(sn,'nhosts')
                    LineConsole.step(['layered queueing network: %d hosts, %d tasks, ' ...
                        '%d entries, %d activities, %d calls'], sn.nhosts, sn.ntasks, ...
                        sn.nentries, sn.nacts, sn.ncalls);
                end
                return
            end
            kind = LineConsole.modelKind(sn);
            LineConsole.step('recognized %s %s: %s, %s, %s', ...
                LineConsole.article(kind), kind, LineConsole.plural(sn.nstations,'station'), ...
                LineConsole.plural(sn.nclasses,'class','classes'), ...
                LineConsole.plural(sn.nchains,'chain'));
            LineConsole.substep('scheduling: %s', LineConsole.schedMix(sn));
            pop = LineConsole.populationLine(sn);
            if ~isempty(pop)
                LineConsole.substep('populations: %s', pop);
            end
            feats = LineConsole.featureList(self);
            if ~isempty(feats)
                LineConsole.substep('features in use: %s', feats);
            end
        end

        function reportMethod(self, options)
            % REPORTMETHOD(SELF, OPTIONS) the method that will actually run
            resolved = options.method;
            if ismethod(self,'resolveMethod')
                resolved = self.resolveMethod(options);
            end
            if strcmp(resolved, options.method)
                LineConsole.step('method ''%s'', tolerance %g, iteration cap %d', ...
                    resolved, options.tol, options.iter_max);
            else
                LineConsole.step('method ''%s'' resolves to ''%s'', tolerance %g, iteration cap %d', ...
                    options.method, resolved, options.tol, options.iter_max);
            end
            if any(strcmp(LineConsole.state().tag, {'SSA','LDES','JMT'})) && ...
                    isfield(options,'samples') && isfinite(options.samples)
                LineConsole.substep('sample budget %g, seed %d', options.samples, options.seed);
            end
        end

        function presolve(self)
            % PRESOLVE(SELF) demands, bottleneck and the elementary bounds
            sn = LineConsole.structOf(self);
            if isempty(sn) || ~LineConsole.has(sn,'nstations') || ...
                    sn.nchains == 0 || isempty(sn.visits) || sn.nstations == 0
                return
            end
            LineConsole.step('computing service demands per chain');
            [Lchain,~,~,~,Nchain,~,~] = sn_get_demands_chain(sn);
            isDelay = LineConsole.delayMask(sn);
            % a Source holds no queue and an infinite server never saturates,
            % so neither is a bottleneck candidate; only the delays are think
            % time, so the two masks stay apart
            notQueueing = isDelay | LineConsole.sourceMask(sn);
            for c = 1:sn.nchains
                L = Lchain(:,c);
                L(~isfinite(L)) = 0;
                Z = sum(L(isDelay));
                Lq = L; Lq(notQueueing) = 0;
                [Dmax, bidx] = max(Lq);
                if Dmax <= GlobalConstants.Zero
                    LineConsole.substep('chain %d has no queueing demand (delay only)', c);
                    continue
                end
                bname = LineConsole.stationName(sn, bidx);
                if isfinite(Nchain(c)) && Nchain(c) > 0
                    LineConsole.substep(['chain %d closed, N = %g, Z = %g: bottleneck %s ' ...
                        'at D = %g, so X <= %g and the knee is at N* = %.3g'], ...
                        c, Nchain(c), Z, bname, Dmax, 1/Dmax, (sum(Lq)+Z)/Dmax);
                else
                    lambda = LineConsole.openRate(sn, c);
                    LineConsole.substep(['chain %d open, lambda = %g: bottleneck %s ' ...
                        'at D = %g, utilization %.4g'], c, lambda, bname, Dmax, lambda*Dmax);
                end
            end
        end

        %% ---------------- the closing narration ----------------

        function closingLines(self)
            % CLOSINGLINES(SELF) the figures of merit and the timing split
            st = LineConsole.state();
            res = LineConsole.resultOf(self);
            if isempty(res)
                LineConsole.ensembleClosing(self, st);
                return
            end
            iter = LineConsole.field(res,'iter',NaN);
            method = LineConsole.field(res,'method','default');
            mtype = line_method_type(st.tag, method);
            if isnan(iter) || iter <= 1
                LineConsole.step('solved by %s (%s)', method, mtype);
            else
                LineConsole.step('solved by %s (%s) in %d iterations', method, mtype, round(iter));
            end
            LineConsole.resultLines(self, res);
            total = toc(st.t0);
            if isnan(st.tsetup)
                LineConsole.step('DONE in %.4f s', total);
            else
                LineConsole.step('DONE in %.4f s (setup %.4f s, analysis %.4f s)', ...
                    total, st.tsetup, max(0,total-st.tsetup));
            end
        end

        function resultLines(self, res)
            % RESULTLINES(SELF, RES) one line per figure of merit
            sn = LineConsole.structOf(self);
            X = LineConsole.field(res,'X',[]);
            C = LineConsole.field(res,'C',[]);
            if ~isempty(X) && ~isempty(C)
                LineConsole.substep('system throughput %s, system response time %s', ...
                    LineConsole.vec(sum(X,1)), LineConsole.vec(sum(C,1)));
            elseif ~isempty(X)
                LineConsole.substep('system throughput %s', LineConsole.vec(sum(X,1)));
            end
            U = LineConsole.field(res,'U',[]);
            if ~isempty(U) && ~isempty(sn) && LineConsole.has(sn,'nstations')
                % a delay station reports jobs in service, not a busy
                % fraction, and a Source has no server at all, so neither is
                % a candidate for the utilization maximum
                Ust = sum(U,2);
                queueing = ~(LineConsole.delayMask(sn) | LineConsole.sourceMask(sn));
                queueing = queueing(1:min(numel(queueing),numel(Ust)));
                if any(queueing)
                    Uq = Ust(1:numel(queueing));
                    Uq(~queueing) = -Inf;
                    [umax, uidx] = max(Uq);
                    LineConsole.substep('busiest queueing station %s at utilization %.4f', ...
                        LineConsole.stationName(sn, uidx), umax);
                end
            end
            Q = LineConsole.field(res,'Q',[]);
            if ~isempty(Q)
                LineConsole.substep('mean jobs in the network %.4f', sum(Q(:)));
            end
        end

        function ensembleClosing(self, st)
            % ENSEMBLECLOSING(SELF, ST) closing lines of an ensemble run
            if isprop(self,'maxitererr') && ~isempty(self.maxitererr)
                err = self.maxitererr(self.maxitererr > 0);
                if isempty(err)
                    err = 0;
                end
                LineConsole.step(['fixed point reached after %d iterations, ' ...
                    'final error %.3e against tolerance %.3e'], ...
                    numel(self.maxitererr), err(end), self.options.iter_tol);
            else
                LineConsole.step('solved');
            end
            LineConsole.step('DONE in %.4f s', toc(st.t0));
        end

        %% ---------------- analyzer hooks ----------------

        function step(fmt, varargin)
            % STEP(FMT, ...) writes one progress line
            if ~LineConsole.writes()
                return
            end
            LineConsole.emitLine('', sprintf(fmt, varargin{:}));
        end

        function substep(fmt, varargin)
            % SUBSTEP(FMT, ...) writes one indented progress line
            if ~LineConsole.writes()
                return
            end
            LineConsole.emitLine('  ', sprintf(fmt, varargin{:}));
        end

        function compileDetail(fmt, varargin)
            % COMPILEDETAIL(FMT, ...) one stage line of a structure compile
            %
            % Silenced inside a PUSHQUIET scope, and inside an open run, where
            % the structures being compiled are those of auxiliary models
            % (SolverLN layers) rather than of the model under study. BEGINRUN
            % lifts the second rule while it compiles its own model.
            st = LineConsole.state();
            if st.quiet > 0 || (st.depth > 0 && ~st.forceDetail)
                return
            end
            LineConsole.substep(fmt, varargin{:});
        end

        function compiling(name)
            % COMPILING(NAME) announces the compilation of a model structure
            %
            % Inside an open run the model being compiled is an auxiliary one
            % (a SolverLN layer), and saying so keeps it apart from the model
            % the user asked about.
            st = LineConsole.state();
            if st.depth > 0 && ~strcmp(name, st.modelName)
                % an ensemble rebuilds the same submodel once per stage or per
                % iteration, so these go through DETAIL and collapse to one
                % line while the name does not change
                LineConsole.detail(sprintf('refreshing the auxiliary submodel ''%s''', name));
            else
                if st.depth == 0 % a compile outside any run opens its own timeline
                    LineConsole.sessionClock(true);
                end
                LineConsole.step('compiling the network structure of model ''%s''', name);
            end
        end

        function detail(text, reset)
            % DETAIL(TEXT) reports a solver's own debug message as a substep
            %
            % LINE_DEBUG routes here while a run narrates. Consecutive repeats
            % are dropped, at most three messages of the same SHAPE (the text
            % with its numbers masked) are reported, and the channel is capped
            % per run, since a message inside a loop would bury the narration.
            %
            % DETAIL('', TRUE) resets the bookkeeping at the start of a run.
            % It is held in persistent variables of this method rather than in
            % the shared state struct: growing a cell array nested inside that
            % struct, thousands of times through a static method, crashed the
            % MATLAB JIT thread outright (segmentation violation) on a long
            % SolverENV run.
            persistent LASTTEXT SHAPES SHAPECOUNT TOTAL
            if isempty(TOTAL)
                LASTTEXT = ''; SHAPES = {}; SHAPECOUNT = []; TOTAL = 0;
            end
            if nargin > 1 && reset
                LASTTEXT = ''; SHAPES = {}; SHAPECOUNT = []; TOTAL = 0;
                return
            end
            if ~LineConsole.ownsLog()
                return
            end
            text = strtrim(text); % routed messages often end in a newline
            if isempty(text) || strcmp(text, LASTTEXT)
                return
            end
            shape = regexprep(text, '[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?', '#');
            hit = find(strcmp(shape, SHAPES), 1);
            if isempty(hit)
                SHAPES{end+1} = shape;
                SHAPECOUNT(end+1) = 1;
            else
                SHAPECOUNT(hit) = SHAPECOUNT(hit) + 1;
                if SHAPECOUNT(hit) > 3
                    return
                end
            end
            LASTTEXT = text;
            TOTAL = TOTAL + 1;
            if TOTAL > 200
                if TOTAL == 201
                    LineConsole.emitLine('  ', 'further solver detail not reported');
                end
                return
            end
            % a message that already reads as a sentence keeps its own wording
            LineConsole.substep('%s', LineConsole.lowerFirst(text));
        end

        function s = lowerFirst(s)
            % S = LOWERFIRST(S) lowercase the first letter of a sentence
            %
            % The routed messages were written as standalone sentences; the
            % console reads as one narration, so they join it in lower case
            % unless they open with an acronym (CTMC, AMVA, LQN, ...).
            if numel(s) >= 2 && ~strcmp(s(1:2), upper(s(1:2)))
                s(1) = lower(s(1));
            end
        end

        function loop(fmt, varargin)
            % LOOP(FMT, ...) announces an iteration loop and resets its budget
            %
            % Each loop of a run gets its own budget of reported iterations,
            % so a solver that restarts a loop (the fluid integration passes)
            % does not exhaust the budget of the next one.
            if ~LineConsole.ownsLog()
                return
            end
            text = sprintf(fmt, varargin{:});
            st = LineConsole.state();
            if strcmp(text, st.lastLoop)
                % the same loop restarted (the fluid integration passes):
                % keep the budget it already spent and say nothing again
                return
            end
            st.lastLoop = text;
            st.shown = 0;
            st.trunc = false;
            LineConsole.state(st);
            LineConsole.emitLine('', text);
        end

        function iter(k, fmt, varargin)
            % ITER(K, FMT, ...) reports iteration K of the current loop
            %
            % Lines are decimated: the first 20 iterations report in full,
            % then every 10th, and the loop stops reporting after 30 lines, so
            % that a long run cannot bury the rest of the narration.
            if ~LineConsole.ownsLog()
                return
            end
            if k > 20 && mod(k,10) ~= 0
                return
            end
            st = LineConsole.state();
            if st.shown >= 30
                if ~st.trunc
                    st.trunc = true;
                    LineConsole.state(st);
                    LineConsole.emitLine('  ', 'further iterations of this loop not reported');
                end
                return
            end
            st.shown = st.shown + 1;
            LineConsole.state(st);
            LineConsole.emitLine('  ', sprintf(fmt, varargin{:}));
        end

        %% ---------------- formatting helpers ----------------

        function t = sessionClock(restart)
            % T = SESSIONCLOCK(RESTART) elapsed seconds on the console clock
            %
            % The clock is restarted by the opening line of each run, so that
            % every run's timeline starts at zero rather than accumulating
            % across the session.
            persistent SESSIONCLOCK
            if (nargin > 0 && restart) || isempty(SESSIONCLOCK)
                SESSIONCLOCK = tic;
            end
            t = toc(SESSIONCLOCK);
        end

        function emitLine(indent, text)
            % EMITLINE(INDENT, TEXT) writes one timestamped line
            %
            % The timestamp is the elapsed time since the opening line of the
            % current run. Each run also reports its own duration in its
            % closing line. A top-level row opens with a capital, an indented
            % substep stays lowercase.
            if isempty(indent) && ~isempty(text)
                text(1) = upper(text(1));
            end
            fprintf(GlobalConstants.StdOut, '[%8.3fs] %s%s\n', LineConsole.sessionClock(), indent, text);
        end

        function s = article(word)
            % S = ARTICLE(WORD) the indefinite article that fits WORD
            if any(lower(word(1)) == 'aeiou')
                s = 'an';
            else
                s = 'a';
            end
        end

        function s = plural(n, singular, plural)
            % S = PLURAL(N, SINGULAR, PLURAL) count with an agreeing noun
            if nargin < 3
                plural = [singular 's'];
            end
            if n == 1
                s = sprintf('%d %s', n, singular);
            else
                s = sprintf('%d %s', n, plural);
            end
        end

        function s = vec(v)
            % S = VEC(V) compact rendering of a numeric row
            v = v(:)';
            if numel(v) == 1
                s = sprintf('%.4f', v);
                return
            end
            parts = cell(1,numel(v));
            for k = 1:numel(v)
                parts{k} = sprintf('%.4f', v(k));
            end
            s = ['[' strjoin(parts,' ') ']'];
        end

        %% ---------------- model introspection ----------------

        function tag = solverTag(self)
            % TAG = SOLVERTAG(SELF) short solver name, e.g. MVA
            tag = strrep(class(self),'Solver','');
            if isempty(tag)
                tag = class(self);
            end
        end

        function model = modelOf(self)
            % MODEL = MODELOF(SELF) the analyzed model, or []
            model = [];
            if isprop(self,'model') && ~isempty(self.model)
                model = self.model;
            end
        end

        function name = modelName(self)
            % NAME = MODELNAME(SELF) the analyzed model's name
            name = '(unnamed)';
            model = LineConsole.modelOf(self);
            if ~isempty(model) && ismethod(model,'getName')
                nm = model.getName();
                if ~isempty(nm)
                    name = char(nm);
                end
            end
        end

        function sn = structOf(self)
            % SN = STRUCTOF(SELF) the model's struct, or [] when it has none
            sn = [];
            model = LineConsole.modelOf(self);
            if ~isempty(model) && ismethod(model,'getStruct')
                sn = model.getStruct();
            end
        end

        function res = resultOf(self)
            % RES = RESULTOF(SELF) the Avg result block, or []
            res = [];
            if isprop(self,'result') && ~isempty(self.result) && ...
                    isfield(self.result,'Avg')
                res = self.result.Avg;
            end
        end

        function tf = has(sn, name)
            % TF = HAS(SN, NAME) field test that works on structs and objects
            if isstruct(sn)
                tf = isfield(sn, name);
            else
                tf = isprop(sn, name);
            end
        end

        function v = field(s, name, dflt)
            % V = FIELD(S, NAME, DFLT) struct field with a default
            if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
                v = s.(name);
            else
                v = dflt;
            end
        end

        function s = versionString()
            % S = VERSIONSTRING() the running LINE version
            v = GlobalConstants.Version;
            if isempty(v)
                s = '';
            else
                s = char(v);
            end
        end

        function s = modelKind(sn)
            % S = MODELKIND(SN) the phrase that names this kind of model
            if any(sn.nodetype == NodeType.Transition) || any(sn.nodetype == NodeType.Place)
                s = 'stochastic Petri net';
                return
            end
            if any(sn.nodetype == NodeType.Cache)
                base = 'caching network';
            elseif any(sn.nodetype == NodeType.Fork)
                base = 'fork-join network';
            else
                base = 'queueing network';
            end
            nopen = sum(isinf(sn.njobs(:)));
            if nopen == 0
                s = ['closed ' base];
            elseif nopen == sn.nclasses
                s = ['open ' base];
            else
                s = ['mixed ' base];
            end
        end

        function s = schedMix(sn)
            % S = SCHEDMIX(SN) counts of each scheduling strategy in use
            scheds = unique(sn.sched(:))';
            parts = {};
            for t = scheds
                n = sum(sn.sched == t);
                parts{end+1} = sprintf('%s x%d', upper(SchedStrategy.toText(t)), n); %#ok<AGROW>
            end
            s = strjoin(parts, ', ');
        end

        function s = populationLine(sn)
            % S = POPULATIONLINE(SN) per-class population or open marker
            parts = {};
            for r = 1:sn.nclasses
                if isinf(sn.njobs(r))
                    parts{end+1} = sprintf('%s open', LineConsole.className(sn,r)); %#ok<AGROW>
                else
                    parts{end+1} = sprintf('%s N=%g', LineConsole.className(sn,r), sn.njobs(r)); %#ok<AGROW>
                end
            end
            s = strjoin(parts, ', ');
        end

        function mask = delayMask(sn)
            % MASK = DELAYMASK(SN) stations that serve without queueing
            mask = false(sn.nstations,1);
            for i = 1:sn.nstations
                mask(i) = sn.sched(i) == SchedStrategy.INF;
            end
        end

        function mask = sourceMask(sn)
            % MASK = SOURCEMASK(SN) stations that inject rather than serve
            mask = false(sn.nstations,1);
            for i = 1:sn.nstations
                mask(i) = sn.sched(i) == SchedStrategy.EXT;
            end
        end

        function lambda = openRate(sn, c)
            % LAMBDA = OPENRATE(SN, C) total arrival rate of an open chain
            lambda = 0;
            inchain = sn.inchain{c};
            for r = inchain(:)'
                for i = 1:sn.nstations
                    if sn.sched(i) == SchedStrategy.EXT && ~isnan(sn.rates(i,r))
                        lambda = lambda + sn.rates(i,r);
                    end
                end
            end
        end

        function s = stationName(sn, i)
            % S = STATIONNAME(SN, I) display name of a station index
            if LineConsole.has(sn,'nodenames') && LineConsole.has(sn,'stationToNode') && ...
                    i >= 1 && i <= numel(sn.stationToNode)
                s = char(sn.nodenames{sn.stationToNode(i)});
            else
                s = sprintf('station %d', i);
            end
        end

        function s = className(sn, r)
            % S = CLASSNAME(SN, R) display name of a class index
            if LineConsole.has(sn,'classnames') && r <= numel(sn.classnames)
                s = char(sn.classnames{r});
            else
                s = sprintf('class %d', r);
            end
        end

        function s = featureList(self)
            % S = FEATURELIST(SELF) the model's used language features
            s = '';
            model = LineConsole.modelOf(self);
            if isempty(model) || ~ismethod(model,'getUsedLangFeatures')
                return
            end
            used = model.getUsedLangFeatures();
            if ~isobject(used) || ~isprop(used,'list')
                return
            end
            names = fieldnames(used.list);
            on = {};
            for k = 1:numel(names)
                if used.list.(names{k})
                    on{end+1} = names{k}; %#ok<AGROW>
                end
            end
            % the structural names (node, class and distribution kinds) only
            % restate the lines printed just above, so they are left out
            structural = {'Source','Sink','Queue','Delay','Router','ClassSwitch', ...
                'OpenClass','ClosedClass','Fork','Join','Cache','Place','Transition', ...
                'Exp','Erlang','HyperExp','Coxian','APH','PH','Det','Immediate','Disabled', ...
                'SchedStrategy_FCFS','SchedStrategy_PS','SchedStrategy_INF','SchedStrategy_EXT', ...
                'RoutingStrategy_PROB','RoutingStrategy_RAND'};
            on = setdiff(on, structural);
            if numel(on) > 8
                s = sprintf('%s, and %d more', strjoin(on(1:8),', '), numel(on)-8);
            else
                s = strjoin(on, ', ');
            end
        end

    end
end
