classdef LineResultRecorder < handle
    % LINERESULTRECORDER Capture result tables with the solver that produced them.
    %
    %   Cross-codebase parity is asserted against one shared golden per example
    %   (goldens/baselines/*.json), keyed by SOLVER NAME. Until 2026-08-19 the
    %   only way to recover that key was to scrape the banner an example printed
    %   above each table -- one regex dialect per codebase, and a value truncated
    %   to whatever the printer showed. The recorder supplies the same
    %   attribution BY CONSTRUCTION: a getter knows which solver it belongs to,
    %   which method that solver resolved, and the full-precision values it is
    %   returning. This is the MATLAB twin of python/line_solver/result_recorder.py.
    %
    %   It is OFF unless asked for, and an ordinary run pays one appdata lookup
    %   per getter and nothing else -- no object is even created:
    %
    %       LineResultRecorder.enable();    % in-process
    %       setenv('LINE_RECORD_RESULTS','1'); LineResultRecorder.autoEnable();
    %
    %   WHAT IS RECORDED. One entry per OUTERMOST getter call. SolverAUTO and the
    %   ensemble solvers call a member solver's getter of the same name, and
    %   recording both would report the member where the example named the
    %   ensemble; the depth guard in ENTER/CAPTURE keeps the outermost only.
    %
    %   WHY THE STATE LIVES IN ROOT APPDATA. Fifteen of the examples in the
    %   corpus run `clear all`, which wipes globals and persistents alike. A
    %   buffer stored in either would come back empty half way through those
    %   runs and read downstream as a solver that produced nothing -- a parity
    %   failure with no defect behind it. Root (handle 0) appdata survives
    %   `clear all` and `close all force`.
    %
    %   See also LINERESULTRECORDER.ENTER, LINERESULTRECORDER.CAPTURE.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant, Hidden)
        % Root appdata key. One recorder per MATLAB session.
        APPDATA_KEY = 'LINEResultRecorder';

        % Solver class name -> the label the shared goldens key it by. These are
        % the names the printed banners carried, which is what the goldens were
        % generated from, so this map is a transcription and not a new scheme.
        SOLVER_LABELS = { ...
            'SolverMVA',  'MVA'; ...
            'SolverNC',   'NC'; ...
            'SolverCTMC', 'CTMC'; ...
            'SolverSSA',  'SSA'; ...
            'SolverFLD',  'FLD'; ...
            'SolverMAM',  'MAM'; ...
            'SolverJMT',  'JMT'; ...
            'SolverLDES', 'LDES'; ...
            'SolverLQNS', 'LQNS'; ...
            'SolverQNS',  'QNS'; ...
            'SolverBA',   'BA'; ...
            'SolverAG',   'AG'; ...
            'SolverAUTO', 'AUTO'; ...
            'SolverLN',   'LN'; ...
            'SolverENV',  'ENV'; ...
            'UQ',         'UQ'};

        % Solvers that RUN ANOTHER SOLVER, and are therefore qualified by it.
        % The member is not decoration: MVA layers and NC layers are different
        % fixed points, and an environment answer scored against one stage's
        % standalone table would report the coupling itself as a defect.
        ENSEMBLE_SOLVERS = {'LN', 'ENV', 'UQ'};

        % The method-family method name FINDSOLVER reports -> the label the goldens key
        % that family by. SolverAUTO(model, 'mva.schmidt') is the documented way
        % to ask for one runnable (family, method) pair -- it is literally the
        % Method column of MODEL.FINDSOLVER -- and a solve asked for that way is
        % recorded under the QUALIFIED label 'MVA:schmidt', which is the golden
        % key for a non-default method.
        %
        % ONLY THE AUTO SPELLING IS QUALIFIED, and that is what keeps the 209
        % existing baselines untouched: SolverMVA(model,'exact') still records as
        % plain 'MVA', because a bare golden key means "this family's default
        % path AS THE EXAMPLE DROVE IT" and cqn_oneline pins 'exact' while keying
        % 'MVA'. Re-labelling direct constructions would move every such golden
        % onto a key no example produces.
        METHOD_FAMILY_LABELS = { ...
            'mva',   'MVA'; ...
            'nc',    'NC'; ...
            'ctmc',  'CTMC'; ...
            'fluid', 'FLD'; ...
            'mam',   'MAM'; ...
            'ba',    'BA'; ...
            'ssa',   'SSA'; ...
            'ldes',  'LDES'; ...
            'jmt',   'JMT'; ...
            'ln',    'LN'; ...
            'env',   'ENV'; ...
            'uq',    'UQ'; ...
            'lqns',  'LQNS'; ...
            'qns',   'QNS'; ...
            'ag',    'AG'};

        % Columns that are never metrics. Anything else numeric in a recorded
        % table is kept, so a solver reporting a metric the goldens do not carry
        % is recorded rather than dropped -- the comparator iterates the GOLDEN,
        % so an extra column costs nothing and a missing one cannot be recovered.
        LABEL_COLUMNS = {'Station', 'Node', 'JobClass', 'Chain', 'Item', ...
                         'Cache', 'Region', 'Class', 'Name', 'NodeType', 'Type'};
    end

    properties (Access = public)
        enabled = false;    % recording on/off
        depth = 0;          % getter nesting depth; only depth 1 records
        records = {};       % cell array of record structs
        notes = {};         % non-table facts a run reported (refusals)
        seq = 0;            % call order, 0-based, as python's recorder numbers
    end

    methods (Access = private)
        function self = LineResultRecorder()
            % Private: use LineResultRecorder.instance().
        end
    end

    methods (Static)
        function r = instance()
            % R = INSTANCE() The session recorder, created on first use.
            r = getappdata(0, LineResultRecorder.APPDATA_KEY);
            if isempty(r) || ~isa(r, 'LineResultRecorder') || ~isvalid(r)
                r = LineResultRecorder();
                setappdata(0, LineResultRecorder.APPDATA_KEY, r);
            end
        end

        function tf = isEnabled()
            % TF = ISENABLED() True when a getter should record.
            %
            % Deliberately does NOT create the recorder: a run that never asked
            % to record must not pay for an object it will never read.
            r = getappdata(0, LineResultRecorder.APPDATA_KEY);
            tf = ~isempty(r) && isa(r, 'LineResultRecorder') && isvalid(r) && r.enabled;
        end

        function enable()
            % ENABLE() Start recording, and drop anything recorded before.
            r = LineResultRecorder.instance();
            r.enabled = true;
            r.reset();
        end

        function disable()
            % DISABLE() Stop recording. The buffer is left intact.
            r = getappdata(0, LineResultRecorder.APPDATA_KEY);
            if ~isempty(r) && isa(r, 'LineResultRecorder') && isvalid(r)
                r.enabled = false;
            end
        end

        function tf = autoEnable()
            % TF = AUTOENABLE() Enable if LINE_RECORD_RESULTS is set to 1.
            tf = strcmp(getenv('LINE_RECORD_RESULTS'), '1');
            if tf
                LineResultRecorder.enable();
            end
        end

        function [depth, guard] = enter()
            % [DEPTH, GUARD] = ENTER() Open a getter's recording scope.
            %
            % DEPTH is 0 when recording is off, which makes CAPTURE a no-op and
            % keeps an ordinary run free of everything below. GUARD is an
            % ONCLEANUP whose destruction pops the depth, so it is popped on
            % EVERY exit from the getter -- including one that raises, which is
            % how the guard survives an example that catches a solver error and
            % carries on.
            %
            % THE GUARD IS RETURNED SEPARATELY, not folded into a struct with
            % the depth, because MATLAB does not promise to destroy an ONCLEANUP
            % held inside a struct or a cell when the enclosing value goes out
            % of scope. A depth left elevated would make every later getter in
            % the same run look like a nested call and record nothing.
            depth = 0;
            guard = [];
            r = getappdata(0, LineResultRecorder.APPDATA_KEY);
            if isempty(r) || ~isa(r, 'LineResultRecorder') || ~isvalid(r) || ~r.enabled
                return
            end
            r.depth = r.depth + 1;
            depth = r.depth;
            guard = onCleanup(@() LineResultRecorder.pop());
        end

        function pop()
            % POP() Close one recording scope. Called by ENTER's onCleanup.
            r = getappdata(0, LineResultRecorder.APPDATA_KEY);
            if ~isempty(r) && isa(r, 'LineResultRecorder') && isvalid(r)
                r.depth = max(0, r.depth - 1);
            end
        end

        function capture(depth, solver, view, tbl)
            % CAPTURE(DEPTH, SOLVER, VIEW, TBL) Record one returned table.
            %
            % DEPTH is ENTER's first output; 0 (recording off) and any nested
            % depth record nothing. VIEW names which table this is.
            if depth ~= 1
                return
            end
            label = LineResultRecorder.solverLabel(solver);
            if isempty(label)
                return
            end
            [rows, labels] = LineResultRecorder.tableRows(tbl, view);
            if isempty(rows)
                return
            end
            r = LineResultRecorder.instance();
            r.append(label, LineResultRecorder.solverMethod(solver), view, ...
                     labels, rows, false);
        end

        function captureScalar(depth, solver, quantity, value)
            % CAPTURESCALAR(DEPTH, SOLVER, QUANTITY, VALUE) Record a derived scalar.
            %
            % Twenty-nine of the goldens hold a quantity no result table carries
            % -- a state probability, a workflow's phase-type moments, a cache
            % hit ratio. Where a library call returns the value, the same
            % construction that records a table records it: the value, the call
            % that produced it, and the solver it belongs to. ParityDerived is
            % what maps a recorded scalar onto the key a particular golden uses,
            % because THAT part is example-specific.
            %
            % The value is flattened to a list of doubles because these getters
            % return a scalar, a pair or an array depending on the solver; the
            % consumer names the element its golden means. DEPTH is ENTER's
            % first output, so a delegating solver records once and not twice.
            if depth ~= 1
                return
            end
            key = LineResultRecorder.solverLabel(solver);
            if isempty(key)
                return
            end
            values = LineResultRecorder.asDoubles(value);
            if isempty(values)
                return
            end
            rows = cell(1, numel(values));
            for i = 1:numel(values)
                rows{i} = struct('Quantity', quantity, ...
                                 'Index', num2str(i - 1), ...
                                 'QLen', values(i));
            end
            r = LineResultRecorder.instance();
            r.append(key, LineResultRecorder.solverMethod(solver), 'scalar', ...
                     {'Quantity', 'Index'}, rows, true);
        end

        function note(text)
            % NOTE(TEXT) Record something that is not a table.
            %
            % Used for a solver's own refusal ('this engine does not implement
            % X'), which a consumer must be able to tell apart from a table that
            % simply never arrived: the first is a fact about the port and is a
            % named skip, the second is a failure.
            if ~LineResultRecorder.isEnabled()
                return
            end
            r = LineResultRecorder.instance();
            r.notes{end+1} = char(text);
        end

        function out = dump()
            % OUT = DUMP() Everything recorded, as a struct ready for JSONENCODE.
            r = LineResultRecorder.instance();
            out = struct('records', {r.records}, 'notes', {r.notes});
        end

        function resetBuffer()
            % RESETBUFFER() Drop the buffer without changing the on/off flag.
            r = LineResultRecorder.instance();
            r.reset();
        end

        function label = solverLabel(solver)
            % LABEL = SOLVERLABEL(SOLVER) The golden's key for this solver.
            %
            % A LAYERED or ENVIRONMENT solve is qualified by the member solver it
            % drove -- 'LN(NC)', 'ENV(FLD)' -- because that is how several
            % goldens spell it and because the two are genuinely different
            % computations. The comparator still reconciles the qualified and
            % bare spellings against the golden's own key; recording the member
            % is what gives it the evidence to do so safely.
            label = '';
            if ~isobject(solver)
                return
            end
            base = LineResultRecorder.labelForClass(class(solver));
            if isempty(base)
                % A subclass of a known solver keys as its parent, the way the
                % `LINE` alias of SolverAUTO does.
                names = LineResultRecorder.SOLVER_LABELS;
                for i = 1:size(names, 1)
                    if isa(solver, names{i, 1})
                        base = names{i, 2};
                        break
                    end
                end
            end
            if isempty(base)
                return
            end
            % AUTO ASKED FOR ONE PAIR IS THAT PAIR. SolverAUTO(model,
            % 'mva.schmidt') delegates to exactly the family and method the method name
            % names, so recording it under 'AUTO' would key it by the meta-solver
            % rather than by what ran -- and 'AUTO' is never goldened (it reports
            % under whichever family it PICKED, so its table would be scored
            % against another family's numbers).
            if strcmp(base, 'AUTO')
                qualified = LineResultRecorder.qualifiedLabel(solver);
                if ~isempty(qualified)
                    label = qualified;
                else
                    label = base;
                end
                return
            end
            if ~any(strcmp(base, LineResultRecorder.ENSEMBLE_SOLVERS))
                label = base;
                return
            end
            % AN ENSEMBLE'S OWN METHOD NAMES IT WHEN IT HAS ONE, because that is
            % the distinction the goldens draw. `lqn_moment3` solves the SAME
            % model with the same NC layers twice, default and 'moment3', and
            % prints the two under `LN Results:` and `LN(moment3) Results:`; its
            % golden holds the first, and they differ by 360x on T1. Labelling
            % both `LN(NC)` would let the second overwrite the first and report
            % the default solve as a 99.7% error.
            method = LineResultRecorder.solverMethod(solver);
            known = LineResultRecorder.SOLVER_LABELS(:, 2);
            if ~isempty(method) && ~strcmp(method, 'default') && ...
                    ~any(strcmpi(method, known))
                label = sprintf('%s(%s)', base, method);
                return
            end
            member = LineResultRecorder.memberSolver(solver);
            if isempty(member)
                label = base;
            else
                label = sprintf('%s(%s)', base, member);
            end
        end

        function member = memberSolver(solver)
            % MEMBER = MEMBERSOLVER(SOLVER) The layer/stage solver an ensemble ran.
            %
            % Getting this wrong is not a spelling difference -- MVA layers and
            % NC layers are different fixed points (lcq_threehosts: cache hit 0.5
            % against 0.48331). SolverLN, SolverENV and UQ all populate
            % `self.solvers` while solving, so the member is read from the object
            % that actually ran rather than guessed from the factory handle.
            member = '';
            if isprop(solver, 'solvers')
                stages = solver.solvers;
                if ~isempty(stages)
                    if ~iscell(stages)
                        stages = {stages};
                    end
                    for i = 1:numel(stages)
                        if isempty(stages{i}) || ~isobject(stages{i})
                            continue
                        end
                        member = LineResultRecorder.labelForClass(class(stages{i}));
                        if isempty(member)
                            names = LineResultRecorder.SOLVER_LABELS;
                            for k = 1:size(names, 1)
                                if isa(stages{i}, names{k, 1})
                                    member = names{k, 2};
                                    break
                                end
                            end
                        end
                        if ~isempty(member)
                            return
                        end
                    end
                end
            end
            % A FACTORY THAT WAS NEVER APPLIED STILL NAMES THE ENGINE. Under
            % lang='java'/'python'/'cpp' the backend solves the whole ensemble
            % and no MATLAB layer is built, so `self.solvers` stays empty and the
            % label fell back to a bare 'LN' -- which reconciles with no golden
            % keyed 'NC' (lqn_twotasks) or 'LN(NC)' (lqn_ofbiz), losing those
            % rows on all three non-MATLAB codebases at once while the MATLAB row
            % passed. The factory is probed on a throwaway network, the same
            % resolution JLINE.lnLayerSolverType, CPPLINE.lnLayerSolver and
            % PYLINE.lnLayerSolverName already perform for the same reason.
            if isprop(solver, 'solverFactory') && ...
                    isa(solver.solverFactory, 'function_handle')
                try
                    probe = solver.solverFactory(CPPLINE.probeNetwork());
                    member = LineResultRecorder.labelForClass(class(probe));
                    if isempty(member)
                        names = LineResultRecorder.SOLVER_LABELS;
                        for k = 1:size(names, 1)
                            if isa(probe, names{k, 1})
                                member = names{k, 2};
                                break
                            end
                        end
                    end
                    if ~isempty(member)
                        return
                    end
                catch
                    % A factory that will not build on the probe names nothing;
                    % the option fallback below is still worth trying.
                    member = '';
                end
            end
            % Nothing solved yet, and no factory either: fall back to the option
            % that names the engine, which is what the LN and ENV dispatchers
            % themselves read.
            if isprop(solver, 'options') && isstruct(solver.options)
                opts = solver.options;
                for name = {'method', 'solver', 'layerSolver', 'stageSolver'}
                    if ~isfield(opts, name{1}) || ~ischar(opts.(name{1}))
                        continue
                    end
                    parts = strsplit(opts.(name{1}), '.');
                    token = upper(strtrim(parts{end}));
                    if any(strcmp(token, {'MVA', 'NC', 'COMOM', 'CTMC', ...
                                          'FLUID', 'FLD', 'SSA', 'LQNS', 'MAM'}))
                        if strcmp(token, 'FLUID')
                            token = 'FLD';
                        end
                        member = token;
                        return
                    end
                end
            end
        end

        function label = qualifiedLabel(solver)
            % LABEL = QUALIFIEDLABEL(SOLVER) 'MVA:schmidt' for a pinned AUTO, '' else.
            %
            % SolverAUTO SPLITS the method name it was given: RESOLVEMETHODTOKEN puts
            % the family in SELECTIONMODE and the submethod in OPTIONS.METHOD, so
            % the pair has to be re-joined here. Reading OPTIONS.METHOD alone
            % would label the solve by a bare submethod ('schmidt'), which names
            % no family and matches no golden key.
            label = '';
            if ~isprop(solver, 'selectionMode')
                return
            end
            family = solver.selectionMode;
            method = LineResultRecorder.solverMethod(solver);
            if isempty(family) || isempty(method) || strcmp(method, 'default')
                return
            end
            names = LineResultRecorder.METHOD_FAMILY_LABELS;
            for i = 1:size(names, 1)
                if strcmpi(family, names{i, 1})
                    label = sprintf('%s:%s', names{i, 2}, method);
                    return
                end
            end
        end

        function method = solverMethod(solver)
            % METHOD = SOLVERMETHOD(SOLVER) The method this solver resolved.
            %
            % `MAM(dec.source)` and `MAM(inap)` differ by 236% on the same model,
            % so the method is not decoration: it is what makes a recorded table
            % attributable to the run the golden was generated from.
            method = 'default';
            if isobject(solver) && isprop(solver, 'options') && ...
                    isstruct(solver.options) && isfield(solver.options, 'method') && ...
                    ischar(solver.options.method) && ~isempty(solver.options.method)
                method = solver.options.method;
            end
        end

        function [rows, labels] = tableRows(tbl, view)
            % [ROWS, LABELS] = TABLEROWS(TBL, VIEW) Rows of a returned table.
            %
            % Values come out at FULL precision: quantizing to the golden's
            % printed precision is the comparator's job, and doing it here would
            % throw away the digits a future full-precision golden needs. A cell
            % that will not convert stays a string, which is right for a label
            % column and harmless for anything else -- the comparator only reads
            % the metrics its golden names.
            rows = {};
            labels = {};
            data = LineResultRecorder.asTable(tbl);
            if isempty(data) || height(data) == 0
                return
            end
            cols = data.Properties.VariableNames;
            labels = cols(ismember(cols, LineResultRecorder.LABEL_COLUMNS));
            if isempty(labels)
                declared = LineResultRecorder.viewLabels(view);
                labels = declared(ismember(declared, cols));
            end
            rows = cell(1, height(data));
            for i = 1:height(data)
                row = struct();
                for c = 1:numel(cols)
                    value = data{i, c};
                    if iscell(value)
                        value = value{1};
                    end
                    if any(strcmp(cols{c}, labels))
                        row.(cols{c}) = LineResultRecorder.asChar(value);
                    elseif isnumeric(value) || islogical(value)
                        row.(cols{c}) = double(value);
                    else
                        row.(cols{c}) = LineResultRecorder.asChar(value);
                    end
                end
                rows{i} = row;
            end
        end

        function data = asTable(tbl)
            % DATA = ASTABLE(TBL) The MATLAB table behind a returned result.
            %
            % Handles the IndexedTable wrapper (which holds it on `.data`), a
            % bare table, and anything else by declining.
            data = [];
            if isa(tbl, 'IndexedTable')
                data = tbl.data;
            elseif istable(tbl)
                data = tbl;
            end
            if ~istable(data)
                data = [];
            end
        end

        function names = viewLabels(view)
            % NAMES = VIEWLABELS(VIEW) The two label columns a view declares.
            %
            % Only consulted for a table whose columns carry none of the names in
            % LABEL_COLUMNS, so an unlisted view still records correctly.
            switch view
                case 'node',       names = {'Node', 'JobClass'};
                case 'chain',      names = {'Station', 'Chain'};
                case 'nodechain',  names = {'Node', 'Chain'};
                case 'sys',        names = {'Chain', 'JobClass'};
                case 'cache',      names = {'Cache', 'Item'};
                case 'item',       names = {'Item', 'JobClass'};
                case 'regionloss', names = {'Region', 'JobClass'};
                case 'region',     names = {'Region', 'JobClass'};
                otherwise,         names = {'Station', 'JobClass'};
            end
        end

        function out = asDoubles(value)
            % OUT = ASDOUBLES(VALUE) Every number in a returned scalar, in order.
            out = [];
            if isempty(value)
                return
            end
            if isnumeric(value) || islogical(value)
                out = double(value(:))';
            elseif iscell(value)
                for i = 1:numel(value)
                    out = [out, LineResultRecorder.asDoubles(value{i})]; %#ok<AGROW>
                end
            end
        end

        function text = asChar(value)
            % TEXT = ASCHAR(VALUE) A label cell as plain char.
            if ischar(value)
                text = value;
            elseif isstring(value) || iscategorical(value)
                text = char(value);
            elseif isnumeric(value) || islogical(value)
                text = num2str(value);
            else
                text = '';
            end
        end
    end

    methods (Access = public)
        function reset(self)
            % RESET(SELF) Drop the buffer and the nesting depth.
            self.records = {};
            self.notes = {};
            self.depth = 0;
            self.seq = 0;
        end

        function append(self, label, method, view, labels, rows, derived)
            % APPEND(SELF, ...) Add one record to the buffer.
            self.records{end+1} = struct( ...
                'solver', label, 'method', method, 'view', view, ...
                'labels', {labels}, 'rows', {rows}, 'seq', self.seq, ...
                'derived', derived);
            self.seq = self.seq + 1;
        end
    end

    methods (Static, Hidden)
        function label = labelForClass(name)
            % LABEL = LABELFORCLASS(NAME) Exact class-name lookup, '' when absent.
            label = '';
            names = LineResultRecorder.SOLVER_LABELS;
            hit = find(strcmp(name, names(:, 1)), 1);
            if ~isempty(hit)
                label = names{hit, 2};
            end
        end
    end
end
