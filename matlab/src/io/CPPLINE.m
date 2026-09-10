classdef CPPLINE
    % CPPLINE  MATLAB-to-C++ bridge (lang='cpp').
    %
    % Static helpers that serialize a LINE model to the canonical model.json
    % interchange (linemodel_save), run the C++ solver binary (line-cli) as a
    % subprocess, and marshal the JSON answer back into MATLAB matrices.
    % Mirrors PYLINE.m (lang='python') and the JLINE.m dispatch (lang='java'),
    % with the transport being subprocess + JSON rather than py.* or the JVM.
    % It is the MATLAB twin of python/line_solver/solvers/cpp_dispatch.py and
    % follows the same rules, flag for flag.
    %
    % WHAT FALLS BACK AND WHAT DOES NOT. lang='cpp' is an assertion about what
    % produced the numbers, so this bridge never substitutes another engine: a
    % missing binary, a construct the C++ analyzer refuses, a non-zero exit or
    % unparseable output are all errors. Silently answering from the MATLAB
    % solver would hide exactly the class of defect cross-language comparison
    % exists to find.
    %
    % LayeredNetwork models reach the C++ layered solver through the .lqnx
    % interchange (writeXML), which is the LQNS file format and covers every
    % construct that reader carries, including fan-in/fan-out replication.
    % A model whose think time sits on a non-reference task cannot be written
    % as .lqnx at all, so it goes over model.json (linemodel_save) instead.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        % The exponent line-cli smooths min(n_i, S_i) with under --method
        % pnorm: FluidOptions::pstar in
        % cpp/include/line/solvers/fluid/solver_fluid.h. The CLI exposes no
        % flag for it, so it is the ONLY exponent this bridge can ask for, and
        % pnormMethod refuses any other rather than answering at this one.
        CLI_PSTAR = 20;
    end

    methods (Static)

        %% ---- binary resolution and transport ----------------------------

        function binary = findLineCli()
            % BINARY = FINDLINECLI()
            % Locate the line-cli executable. Search order, most explicit
            % first: the LINE_CLI_BINARY environment variable, the checkout's
            % common/ directory, the in-tree cpp/build directories, then PATH.
            % Errors rather than returning '' so a missing binary cannot be
            % mistaken for a solver refusal further down.
            binary = '';
            envBin = getenv('LINE_CLI_BINARY');
            if ~isempty(envBin)
                if exist(envBin, 'file') == 2
                    binary = envBin;
                    return
                end
                line_error(mfilename, sprintf(['LINE_CLI_BINARY is set to ''%s'', ' ...
                    'which is not an existing file.'], envBin));
            end

            exeName = 'line-cli';
            if ispc
                exeName = 'line-cli.exe';
            end
            root = lineRootFolder();     % .../matlab
            repo = fileparts(root);      % checkout root
            cands = {fullfile(repo, 'common', exeName), ...
                fullfile(repo, 'cpp', 'build', exeName), ...
                fullfile(repo, 'cpp', 'build-gmp', exeName)};
            for k = 1:numel(cands)
                if exist(cands{k}, 'file') == 2
                    binary = cands{k};
                    return
                end
            end

            % PATH lookup, last because it is the least attributable.
            if ispc
                [st, out] = system(['where ' exeName]);
            else
                [st, out] = system(['command -v ' exeName]);
            end
            if st == 0
                out = strtrim(strtok(out, sprintf('\n')));
                if ~isempty(out) && exist(out, 'file') == 2
                    binary = out;
                    return
                end
            end

            line_error(mfilename, ['the C++ solver binary ''line-cli'' was not found. ' ...
                'Set LINE_CLI_BINARY to its path, put it on PATH, or build it with ' ...
                '''cpp/make.sh''.']);
        end

        function results = runLineCli(binary, args)
            % RESULTS = RUNLINECLI(BINARY, ARGS)
            % Run one line-cli invocation (ARGS a cellstr of argv entries after
            % the binary) and return the jsondecode'd object. This is the seam
            % an in-process binding would replace: nothing above it knows a
            % process was involved.
            % MATLAB's own libstdc++ shadows the host one on LD_LIBRARY_PATH,
            % so a binary built against this toolchain fails to load; see
            % line_native_env.
            cmd = [line_native_env() sprintf('"%s"', binary)];
            for k = 1:numel(args)
                cmd = [cmd ' ' CPPLINE.shellQuote(args{k})]; %#ok<AGROW>
            end
            errFile = [tempname '.err'];
            cleaner = onCleanup(@() CPPLINE.tryDelete(errFile)); %#ok<NASGU>
            [status, out] = system([cmd ' 2> ' CPPLINE.shellQuote(errFile)]);
            errText = '';
            if exist(errFile, 'file') == 2
                errText = strtrim(fileread(errFile));
            end
            if status ~= 0
                % A STALE BINARY IS THE FIRST THING TO RULE OUT: -o json on the
                % model-solving path postdates this bridge's need for it, and a
                % binary predating it refuses with this exact phrase.
                if contains(errText, 'prints only -o readable')
                    line_error(mfilename, sprintf(['''%s'' is a STALE line-cli: it predates ' ...
                        'the -o json output path this bridge needs. Rebuild it (cpp/make.sh) ' ...
                        'or point LINE_CLI_BINARY at a fresh one.'], binary));
                end
                % exit 2 is line-cli's refusal channel, with the construct named
                % in the message. It is relayed verbatim because that message is
                % the port's contract with the caller.
                line_error(mfilename, sprintf('line-cli exited with code %d.\nstderr:\n%s\nstdout:\n%s', ...
                    status, errText, strtrim(out)));
            end
            if ~isempty(errText)
                line_warning(mfilename, 'line-cli: %s\n', errText);
            end
            results = CPPLINE.extractJsonObject(out, binary);
        end

        function obj = extractJsonObject(text, binary)
            % OBJ = EXTRACTJSONOBJECT(TEXT, BINARY)
            % Decode the first complete JSON object in TEXT, ignoring the
            % solver banner line-cli prints ahead of it.
            if nargin < 2
                binary = 'line-cli';
            end
            start = strfind(text, '{');
            if isempty(start)
                % THE SILENT STALE-BINARY SIGNATURE: -s ssa and -s fluid used to
                % accept -o json and print the readable table anyway, so an old
                % build exits 0 with a brace-free table.
                if contains(text, 'Station') && contains(text, 'JobClass')
                    line_error(mfilename, sprintf(['''%s'' answered -o json with the READABLE ' ...
                        'table, which means it predates the shared JSON emitter for the ' ...
                        '-s ssa and -s fluid arms. Rebuild it (cpp/make.sh).'], binary));
                end
                line_error(mfilename, sprintf('no JSON object found in line-cli output:\n%s', strtrim(text)));
            end
            jsonText = text(start(1):end);
            stop = CPPLINE.matchingBrace(jsonText);
            if stop > 0
                jsonText = jsonText(1:stop);
            end
            try
                obj = jsondecode(jsonText);
            catch ME
                line_error(mfilename, sprintf('could not parse line-cli JSON output: %s\n%s', ...
                    ME.message, strtrim(text)));
            end
        end

        function stop = matchingBrace(s)
            % STOP = MATCHINGBRACE(S)
            % Index of the '}' closing the '{' at S(1), scanning outside string
            % literals so a brace inside a name cannot end the object early.
            depth = 0;
            inStr = false;
            escaped = false;
            for k = 1:numel(s)
                c = s(k);
                if inStr
                    if escaped
                        escaped = false;
                    elseif c == '\'
                        escaped = true;
                    elseif c == '"'
                        inStr = false;
                    end
                    continue
                end
                switch c
                    case '"'
                        inStr = true;
                    case '{'
                        depth = depth + 1;
                    case '}'
                        depth = depth - 1;
                        if depth == 0
                            stop = k;
                            return
                        end
                end
            end
            stop = 0;
        end

        function q = shellQuote(tok)
            % Q = SHELLQUOTE(TOK)  Quote one argv entry for the platform shell.
            tok = char(tok);
            if ispc
                q = ['"' strrep(tok, '"', '""') '"'];
            else
                q = ['''' strrep(tok, '''', '''\''''') ''''];
            end
        end

        function tryDelete(f)
            % TRYDELETE  Best-effort temp-file removal.
            if exist(f, 'file') == 2
                delete(f);
            end
        end

        %% ---- solver method name and refusals ----------------------------------

        function token = solverToken(solverName)
            % TOKEN = SOLVERTOKEN(SOLVERNAME)
            % The line-cli -s method name for a LINE solver name. The C++
            % model-solving path's token set is NOT the JAR's: it has 'ag' and
            % 'ba', which the JAR CLI lacks. It does carry 'jmt', which drives
            % the same JMT.jar the MATLAB wrapper drives, and lacks lqns and
            % qns's layered binaries.
            %
            % THE TABLE IS THE AUTHORITY FOR WHAT lang='cpp' ACCEPTS, so a token
            % missing here makes a CLI arm unreachable however complete the port
            % is -- `ag` was: `line-cli -s ag` served -a avg, its four views and
            % -a cdf, and native Python's `_CPP_SOLVER_TOKENS` already carried
            % it, while this switch refused SolverAG outright and every `ag_*`
            % [M2C] parity row read as a named gap in the C++ solver.
            name = upper(char(solverName));
            if strncmp(name, 'SOLVER', 6)
                name = name(7:end);
            end
            switch name
                case 'MVA',   token = 'mva';
                case 'NC',    token = 'nc';
                case 'CTMC',  token = 'ctmc';
                case 'MAM',   token = 'mam';
                case {'FLD','FLUID'}, token = 'fluid';
                case 'SSA',   token = 'ssa';
                case 'BA',    token = 'ba';
                case 'AG',    token = 'ag';
                case 'JMT',   token = 'jmt';
                case 'AUTO',  token = 'auto';
                otherwise
                    line_error(mfilename, sprintf(['lang=''cpp'' is not available for solver ' ...
                        '''%s'': the C++ model-solving path serves ag, auto, ba, ctmc, fluid, ' ...
                        'jmt, mam, mva, nc and ssa. LDES is a JSON subprocess client in every ' ...
                        'codebase and keeps its own path, LQNS and QNS wrap layered binaries ' ...
                        'the port does not carry, and a LayeredNetwork goes through the layered ' ...
                        'path instead. Use lang=''matlab'', lang=''java'' or lang=''python''.'], solverName));
            end
        end

        function cppUnsupported(solverName, method, reason)
            % CPPUNSUPPORTED(SOLVERNAME, METHOD, REASON)
            % Refuse a getter this bridge cannot serve, naming the getter and
            % the reason. THE REASON IS THE POINT: what is unreachable is
            % unreachable for a specific, per-getter cause, and collapsing those
            % into one message would tell a caller to stop asking for something
            % the port may be one arm away from answering.
            name = char(solverName);
            if strncmpi(name, 'Solver', 6)
                name = name(7:end);
            end
            line_error(mfilename, sprintf(['lang=''cpp'' cannot serve Solver%s.%s: %s. ' ...
                'Use lang=''matlab'', lang=''java'' or lang=''python'' for it.'], name, method, reason));
        end

        %% ---- model-solving path -----------------------------------------

        function results = solveViaCpp(solverName, model, options, analysis, flags)
            % RESULTS = SOLVEVIACPP(SOLVERNAME, MODEL, OPTIONS, ANALYSIS, FLAGS)
            % Serialize MODEL to model.json, run line-cli on it and return the
            % parsed answer.
            %
            % Only the knobs line-cli actually honours are forwarded, and only
            % when the caller overrode the MATLAB default: an untouched option
            % lets the port apply its own SolverOptions default rather than
            % having this bridge's default imposed on it. KNOBS ARE ALSO GATED
            % BY METHOD NAME, because line-cli REFUSES an option the chosen solver
            % does not have (--tol on -s ba, --samples outside ssa, --cutoff
            % outside ctmc), so a solver whose default merely differs from the
            % generic one would otherwise be refused for an option nobody set.
            %
            % FLAGS is an optional cellstr of the analysis's OWN argv entries
            % (--tspan for a transient query, --node for a per-node one,
            % --notation for the ODE export). They are arguments of the CALL and
            % not settings of the solver, so the caller that knows the getter
            % assembles them.
            if nargin < 4 || isempty(analysis)
                analysis = 'avg';
            end
            if nargin < 5 || isempty(flags)
                flags = {};
            end
            binary = CPPLINE.findLineCli();
            token = CPPLINE.solverToken(solverName);

            tmpDir = tempname;
            if ~mkdir(tmpDir)
                line_error(mfilename, sprintf('could not create the temporary directory ''%s''.', tmpDir));
            end
            cleaner = onCleanup(@() CPPLINE.rmdirQuiet(tmpDir)); %#ok<NASGU>
            modelPath = fullfile(tmpDir, 'model.json');
            linemodel_save(model, modelPath);
            if exist(modelPath, 'file') ~= 2
                line_error(mfilename, 'linemodel_save did not write a model.json for the C++ solver.');
            end

            args = {'-f', modelPath, '-i', 'json', '-s', token, '-a', char(analysis), '-o', 'json'};
            args = [args, CPPLINE.optionFlags(options, token)];
            for k = 1:numel(flags)
                args{end+1} = char(flags{k}); %#ok<AGROW>
            end
            results = CPPLINE.runLineCli(binary, args);
        end

        function args = optionFlags(options, token)
            % ARGS = OPTIONFLAGS(OPTIONS, TOKEN)
            % The subset of SolverOptions line-cli honours for METHOD NAME, forwarded
            % only where the caller overrode the MATLAB default.
            args = {};
            if isempty(options) || ~isstruct(options)
                return
            end
            % --arith is the one flag with no lang='java' counterpart: it is
            % what the C++ port exists for. Absent, the port runs in double,
            % which is the only setting comparable with the other backends.
            if isfield(options, 'arith') && ~isempty(options.arith) && ...
                    ~strcmpi(char(options.arith), 'double')
                args = [args, {'--arith', char(options.arith)}];
            end
            % NOTHING IS RESOLVED HERE FOR THE FLUID SOLVER: line-cli resolves
            % 'default' itself (fluid_resolve_default_method), and forwarding a
            % name computed on this side would override that resolution with one
            % taken from a different model struct.
            if isfield(options, 'method') && ~isempty(options.method) && ...
                    ~strcmp(char(options.method), 'default')
                args = [args, {'--method', char(options.method)}];
            end
            % ba/ssa/ctmc/jmt have no convergence loop to control and line-cli refuses these flags on them; mam carries tol and iter_max but not iter_tol.
            if ~any(strcmp(token, {'ba', 'ssa', 'ctmc', 'jmt'}))
                tol = CPPLINE.overridden(options, 'tol', 1e-4);
                if ~isempty(tol)
                    args = [args, {'--tol', CPPLINE.num2arg(tol)}];
                end
                iterMax = CPPLINE.overridden(options, 'iter_max', 1000);
                if ~isempty(iterMax) && iterMax > 0
                    args = [args, {'--iter_max', sprintf('%d', round(iterMax))}];
                end
                if ~any(strcmp(token, {'mam', 'ag'}))
                    iterTol = CPPLINE.overridden(options, 'iter_tol', 1e-4);
                    if ~isempty(iterTol)
                        args = [args, {'--iter_tol', CPPLINE.num2arg(iterTol)}];
                    end
                end
            end
            % SolverAG's own knob. maxStates is a TRUNCATION LEVEL and therefore
            % part of the answer rather than a budget, so it travels: a run
            % truncated at the port's 100 when the caller asked for 500 is a
            % different number reported as theirs. The execution backends do
            % NOT travel -- line-cli exposes no flag for them and its AgOptions
            % always sweeps serially -- so a non-serial one is refused rather
            % than dropped, which is the same silent-acceptance rule the CLI
            % applies from its side.
            if strcmp(token, 'ag') && isfield(options, 'config') && isstruct(options.config)
                maxStates = [];
                if isfield(options.config, 'maxStates')
                    maxStates = CPPLINE.overridden(options.config, 'maxStates', 100);
                end
                if ~isempty(maxStates) && maxStates > 0
                    args = [args, {'--max-states', sprintf('%d', round(maxStates))}];
                end
                if isfield(options.config, 'exec') && ~isempty(options.config.exec) && ...
                        ~strcmpi(char(options.config.exec), 'serial')
                    line_error(mfilename, sprintf(['lang=''cpp'' cannot serve SolverAG with ' ...
                        'options.config.exec=''%s'': line-cli exposes no execution-backend flag, ' ...
                        'so the port always sweeps the agents serially. The backends produce the ' ...
                        'SAME iterates, so this refuses only the placement of the work, not an ' ...
                        'answer: use lang=''matlab'', lang=''java'' or lang=''python'' when the ' ...
                        'pool or the remote workers are the point.'], char(options.config.exec)));
                end
            end
            % A sample path has a length and a stream; JSIM has both, and dropping them would run at the port's defaults rather than the caller's seed.
            if any(strcmp(token, {'ssa', 'jmt'}))
                if isfield(options, 'samples') && ~isempty(options.samples) && options.samples > 0
                    args = [args, {'--samples', sprintf('%d', round(options.samples))}];
                end
                % options.seed defaults to a positive value in MATLAB, and
                % line-cli refuses 0 outright, so only a positive one is sent.
                if isfield(options, 'seed') && ~isempty(options.seed) && options.seed > 0
                    args = [args, {'--seed', sprintf('%d', round(options.seed))}];
                end
            end
            % The cutoff bounds an open population inside a state space; it is a
            % SolverCTMC option and reaches every CTMC analysis.
            if strcmp(token, 'ctmc')
                if isfield(options, 'cutoff') && ~isempty(options.cutoff)
                    c = max(options.cutoff(:));
                    if isfinite(c) && c > 0
                        args = [args, {'--cutoff', CPPLINE.num2arg(c)}];
                    end
                end
            end
            % The QRF reduction bounds of SolverBA read their blocking tables
            % and their load-dependent scaling from options.config; without
            % them the port refuses qrf.bas and qrf.rsrd outright, exactly as
            % sn_to_qrf_params does here.
            if strcmp(token, 'ba')
                lvl = CPPLINE.overridden(options, 'level', 2);
                if ~isempty(lvl) && lvl > 0
                    args = [args, {'--level', sprintf('%d', round(lvl))}];
                end
                if isfield(options, 'config') && isstruct(options.config)
                    if isfield(options.config, 'qrf_params') && ~isempty(options.config.qrf_params)
                        args = [args, {'--qrf-params', CPPLINE.qrfParamsJson(options.config.qrf_params)}];
                    end
                    if isfield(options.config, 'qrf_alpha') && ~isempty(options.config.qrf_alpha)
                        args = [args, {'--qrf-alpha', jsonencode(CPPLINE.rowCells(options.config.qrf_alpha))}];
                    end
                end
            end
            % AMVA multiserver rule: a rule that does not reach the solver
            % leaves it answering under 'default' while the caller believes
            % their own choice was applied. jmt simulates the multiserver directly and line-cli refuses the flag on it.
            if ~strcmp(token, 'jmt') && isfield(options, 'config') && isstruct(options.config) && ...
                    isfield(options.config, 'multiserver') && ...
                    ~isempty(options.config.multiserver) && ...
                    ~strcmp(char(options.config.multiserver), 'default')
                args = [args, {'--multiserver', char(options.config.multiserver)}];
            end
        end

        function c = rowCells(A)
            % C = ROWCELLS(A)
            % A numeric matrix as a cell of its rows, so that jsonencode always
            % emits an array of arrays. A bare one-row matrix encodes FLAT,
            % which the CLI reads as a vector and not as a one-row table.
            A = double(A);
            c = cell(size(A, 1), 1);
            for i = 1:size(A, 1)
                c{i} = A(i, :);
            end
        end

        function txt = qrfParamsJson(qp)
            % TXT = QRFPARAMSJSON(QP)
            % options.config.qrf_params as the --qrf-params document. ZM is not
            % sent: the CLI derives it from ZZ, as every codebase now does.
            if ~isstruct(qp)
                line_error(mfilename, 'options.config.qrf_params must be a struct.');
            end
            required = {'f', 'MR', 'BB', 'MM', 'MM1', 'ZZ'};
            for i = 1:numel(required)
                if ~isfield(qp, required{i})
                    line_error(mfilename, sprintf(['options.config.qrf_params must contain ' ...
                        'field ''%s'' for the QRF blocking bounds.'], required{i}));
                end
            end
            doc = struct();
            doc.f = double(qp.f);
            doc.MR = double(qp.MR);
            doc.BB = CPPLINE.rowCells(qp.BB);
            doc.MM = CPPLINE.rowCells(qp.MM);
            doc.MM1 = CPPLINE.rowCells(qp.MM1);
            doc.ZZ = double(qp.ZZ(:)');
            if isfield(qp, 'F') && ~isempty(qp.F)
                doc.F = double(qp.F(:)');
            end
            txt = jsonencode(doc);
        end

        function v = overridden(options, field, defaultValue)
            % V = OVERRIDDEN(OPTIONS, FIELD, DEFAULTVALUE)
            % The option value when the caller moved it off the MATLAB default,
            % [] otherwise.
            v = [];
            if ~isfield(options, field) || isempty(options.(field))
                return
            end
            val = options.(field);
            if ~isnumeric(val) || ~isscalar(val) || ~isfinite(val)
                return
            end
            if abs(val - defaultValue) <= 1e-15 * max(1, abs(defaultValue))
                return
            end
            v = val;
        end

        function s = num2arg(x)
            % S = NUM2ARG(X)  Round-trip-exact decimal form of a scalar knob.
            s = sprintf('%.17g', double(x));
        end

        function rmdirQuiet(d)
            % RMDIRQUIET  Best-effort recursive temp-directory removal.
            if exist(d, 'dir') == 7
                try %#ok<TRYNC>
                    rmdir(d, 's');
                end
            end
        end

        %% ---- average table ------------------------------------------------

        function [QN, UN, RN, TN, AN, WN, runtime] = getAvg(solverName, model, options)
            % [QN,UN,RN,TN,AN,WN,RUNTIME] = GETAVG(SOLVERNAME, MODEL, OPTIONS)
            % Run the C++ avg analysis and reduce it to (nstations x nclasses)
            % matrices in the MATLAB struct's own indexing.
            Tstart = tic;
            results = CPPLINE.solveViaCpp(solverName, model, options, 'avg');
            sn = model.getStruct();
            [QN, UN, RN, TN, AN, WN] = CPPLINE.avgMatrices(results, sn);
            CPPLINE.restoreCacheResults(results, model);
            runtime = toc(Tstart);
        end

        function restoreCacheResults(results, model)
            % RESTORECACHERESULTS(RESULTS, MODEL)
            % Write the C++ per-Cache results back onto the MATLAB Cache nodes.
            %
            % A cache's hit, miss and delayed-hit fractions and its retrieval
            % latency are SOLVER RESULTS, not model state: getAvgCacheTable reads
            % them off the node. Solving in the C++ engine used to leave the
            % node untouched, so the table reported whatever the PREVIOUS solver
            % had written -- on retrieval_simple the MVA and NC tables were both
            % the LDES table, identical to five digits, and nothing in those two
            % rows was being measured at all.
            %
            % ABSENT MUST CLEAR, NOT KEEP. Every field is written on every
            % solve, [] included, mirroring the lang='java' path
            % (SolverMVA/runAnalyzer.m:68-75). Writing only the fields the
            % engine filled is what let one solver's delayed-hit fraction
            % survive into another solver's row.
            % The block rides INSIDE the "avg" payload, not beside it: line-cli
            % merges `extra` into the analysis object, the way it already does
            % for ListCost (line_cli.cpp:448). Reading it at envelope level
            % finds nothing and silently restores no cache result at all.
            if ~isa(model, 'Network')
                return
            end
            entries = {};
            if isstruct(results) && isfield(results, 'avg')
                table = results.avg;
                if isstruct(table) && isfield(table, 'Cache')
                    entries = table.Cache;
                end
            end
            % NO BLOCK MEANS THE SOLVER MEASURED NOTHING, WHICH MUST CLEAR. This
            % used to return here, and a solver with no cache branch then left
            % the node holding the previous solver's split -- which is how the
            % MVA and NC rows printed the LDES table. Falling through with an
            % empty list writes [] to every field, so getAvgCacheTable reports
            % NaN and the row disappears, as it does natively.
            if isstruct(entries)
                entries = num2cell(entries);
            end
            if ~iscell(entries)
                entries = {};
            end
            nodes = model.getNodes;
            changed = false;
            % Driven by the MODEL's Cache nodes, not by the payload's entries, so
            % that a cache the solver said nothing about is CLEARED rather than
            % left alone.
            for ind = 1:numel(nodes)
                if ~isa(nodes{ind}, 'Cache')
                    continue
                end
                % MATCHED BY NAME. `e.node` indexes the CLI's OWN node order,
                % which is not this model's: on retrieval_simple the Cache is 2
                % here and 4 there, so an index match landed on the Sink and
                % restored nothing at all.
                e = [];
                for k = 1:numel(entries)
                    c = entries{k};
                    if isstruct(c) && isfield(c, 'name') && ...
                            strcmp(char(c.name), nodes{ind}.getName())
                        e = c;
                        break
                    end
                end
                hit = []; miss = []; dhit = []; residt = [];
                dhq = []; dhqf = [];
                if ~isempty(e)
                    hit = CPPLINE.cacheField(e, 'HitProb');
                    miss = CPPLINE.cacheField(e, 'MissProb');
                    if isempty(miss) && ~isempty(hit)
                        miss = 1 - hit;
                    end
                    dhit = CPPLINE.cacheField(e, 'DelayedHitProb');
                    residt = CPPLINE.cacheField(e, 'ResidT');
                    % Per-item delayed-hit queue lengths, the two DelayedHitQLen
                    % columns of getAvgItemTable. Only the exact chain computes
                    % them, so they are absent from every other C++ solve and
                    % must then CLEAR, like every field above.
                    dhq = CPPLINE.cacheField(e, 'DelayedHitQLen');
                    dhqf = CPPLINE.cacheField(e, 'DelayedHitQLenFull');
                end
                changed = changed || ~isempty(hit) || ~isempty(nodes{ind}.getHitRatio());
                nodes{ind}.setResultHitProb(hit);
                nodes{ind}.setResultMissProb(miss);
                nodes{ind}.setResultDelayedHitProb(dhit);
                nodes{ind}.setResultResidT(residt);
                nodes{ind}.setResultDelayedHitQLen(dhq, dhqf);
            end
            if changed
                % The visits carrying the hit and miss classes are derived from
                % the split, so the struct has to be refreshed hard or a Router
                % downstream of the cache keeps link()'s uniform guess. Same
                % reason as PYLINE.restoreCacheResults. A CLEAR refreshes too:
                % leaving the derived visits behind is the same staleness in the
                % struct that the clear just removed from the node.
                model.refreshStruct(true);
            end
        end

        function v = cacheField(e, name)
            % V = CACHEFIELD(E, NAME)
            % A per-class cache result as a row vector, [] when the solver
            % computed none. The empty is the point: it CLEARS a stale value.
            v = [];
            if isfield(e, name)
                v = CPPLINE.jsonNumericList(e.(name));
                if ~isempty(v)
                    v = v(:)';
                end
            end
        end

        function [QN, UN, RN, TN, AN, WN] = avgMatrices(results, sn)
            % [QN,UN,RN,TN,AN,WN] = AVGMATRICES(RESULTS, SN)
            % Reduce a line-cli avg payload to station x class metric matrices.
            %
            % ROWS ARE KEYED BY NAME, never by position: the CLI emits the
            % stations in the order the model.json declares them, which is the
            % NODE order, while the matrices are indexed by STATION. Reading
            % them positionally would put a station's numbers on another
            % station's row on any model with a non-station node before a
            % station (a Source, a Router, a ClassSwitch).
            if ~isstruct(results) || ~isfield(results, 'avg')
                line_error(mfilename, sprintf(['line-cli did not answer -a avg (keys: %s).'], ...
                    strjoin(fieldnames(results)', ', ')));
            end
            table = results.avg;
            if ~isstruct(table) || ~isfield(table, 'QLen')
                line_error(mfilename, 'line-cli did not return a structured AvgTable.');
            end

            M = sn.nstations;
            K = sn.nclasses;
            QN = zeros(M, K); UN = zeros(M, K); RN = zeros(M, K);
            TN = zeros(M, K); AN = zeros(M, K); WN = zeros(M, K);

            % node name -> station index, from the struct's own hashing maps
            nodeToStation = sn.nodeToStation(:);
            stationOfName = configureDictionary('string', 'double');
            for nd = 1:sn.nnodes
                st = nodeToStation(nd);
                if st >= 1
                    stationOfName(string(sn.nodenames(nd))) = st;
                end
            end
            classOfName = configureDictionary('string', 'double');
            for r = 1:K
                classOfName(string(sn.classnames(r))) = r;
            end

            stations = CPPLINE.jsonStringList(table.Station);
            jobclasses = CPPLINE.jsonStringList(table.JobClass);
            qlen = CPPLINE.jsonNumericList(table.QLen);
            util = CPPLINE.jsonNumericList(table.Util);
            respt = CPPLINE.jsonNumericList(table.RespT);
            residt = CPPLINE.jsonNumericList(table.ResidT);
            arvr = CPPLINE.jsonNumericList(table.ArvR);
            tput = CPPLINE.jsonNumericList(table.Tput);

            for i = 1:numel(stations)
                sname = string(stations{i});
                cname = string(jobclasses{i});
                if ~isKey(stationOfName, sname) || ~isKey(classOfName, cname)
                    continue
                end
                st = stationOfName(sname);
                cl = classOfName(cname);
                QN(st, cl) = qlen(i);
                UN(st, cl) = util(i);
                RN(st, cl) = respt(i);
                WN(st, cl) = residt(i);
                AN(st, cl) = arvr(i);
                TN(st, cl) = tput(i);
            end

            % FINITE-CAPACITY REGIONS RIDE PAST THE LAST STATION. getAvgNode
            % reads self.result.Avg rows nstations+1 .. nstations+nregions for
            % the FCR pseudo-nodes and getAvgNodeTable then drops any node whose
            % metrics are all zero, so a bridge that returned only the station
            % rows deleted the FCR row from the node table outright
            % (fcr_mm1waitq[M2C], 'row FCR1 missing'). line-cli carries them in
            % the avg payload's FCR block, in the region order the struct
            % declares.
            F = 0;
            if isfield(sn, 'nregions') && ~isempty(sn.nregions)
                F = double(sn.nregions);
            end
            if F > 0 && isfield(table, 'FCR') && isstruct(table.FCR) ...
                    && isfield(table.FCR, 'QLen')
                fq = CPPLINE.jsonNumericList(table.FCR.QLen);
                fu = CPPLINE.jsonNumericList(table.FCR.Util);
                fr = CPPLINE.jsonNumericList(table.FCR.RespT);
                fw = CPPLINE.jsonNumericList(table.FCR.ResidT);
                fa = CPPLINE.jsonNumericList(table.FCR.ArvR);
                ft = CPPLINE.jsonNumericList(table.FCR.Tput);
                QN = [QN; zeros(F, K)]; UN = [UN; zeros(F, K)];
                RN = [RN; zeros(F, K)]; WN = [WN; zeros(F, K)];
                AN = [AN; zeros(F, K)]; TN = [TN; zeros(F, K)];
                for f = 1:F
                    for cl = 1:K
                        j = (f - 1) * K + cl;   % row-major, as the CLI writes it
                        if j > numel(fq)
                            continue
                        end
                        QN(M + f, cl) = fq(j);
                        UN(M + f, cl) = fu(j);
                        RN(M + f, cl) = fr(j);
                        WN(M + f, cl) = fw(j);
                        AN(M + f, cl) = fa(j);
                        TN(M + f, cl) = ft(j);
                    end
                end
            end
        end

        function c = jsonStringList(v)
            % C = JSONSTRINGLIST(V)
            % Normalize a jsondecode'd JSON array of strings to a cellstr. A
            % one-element array decodes to a char row and an all-equal-length
            % array may decode to a char matrix, so neither shape may be assumed.
            if iscell(v)
                c = cellfun(@char, v, 'UniformOutput', false);
                c = c(:)';
            elseif isstring(v)
                c = cellstr(v(:)');
            elseif ischar(v)
                if size(v, 1) > 1
                    c = cellstr(v)';
                else
                    c = {v};
                end
            else
                c = {};
            end
        end

        function x = jsonNumericList(v)
            % X = JSONNUMERICLIST(V)
            % Normalize a jsondecode'd JSON array of numbers to a row vector,
            % mapping a JSON null (which arrives as an empty cell entry) to 0.
            % An INFINITY crosses as the string 'Infinity'/'-Infinity', since
            % JSON has no infinite literal; read as a name and defaulted to 0 it
            % turns a queue at its stability boundary into an idle one.
            if iscell(v)
                x = zeros(1, numel(v));
                for k = 1:numel(v)
                    if ischar(v{k}) || isstring(v{k})
                        x(k) = CPPLINE.jsonNumericScalar(v{k});
                    elseif isempty(v{k}) || ~isnumeric(v{k})
                        x(k) = 0;
                    else
                        x(k) = double(v{k}(1));
                    end
                end
            else
                x = double(v(:)');
            end
        end

        function v = jsonNumericScalar(x)
            % V = JSONNUMERICSCALAR(X)
            % One wire value as a double: the non-finite string spellings both
            % CLIs emit, or str2double for anything else numeric-looking.
            %
            % A NUMBER IS RETURNED AS ITSELF, and reaching char() with one was a
            % silent zero: jsondecode gives a double for an ordinary JSON
            % number, char(0.2231) is the character at CODE POINT 0, str2double
            % of that empty string is NaN, and the NaN branch below floored it
            % to 0. Every scalar read through here -- ProbSys, ProbSysAggr,
            % logNormConstAggr -- therefore came back 0 whatever the engine
            % said, and 0 is a legal probability, so nothing complained.
            if isnumeric(x) || islogical(x)
                if isempty(x)
                    v = 0;
                else
                    v = double(x(1));
                end
                return
            end
            s = strtrim(char(x));
            switch s
                case {'Infinity','inf','Inf'}
                    v = Inf;
                case {'-Infinity','-inf','-Inf'}
                    v = -Inf;
                case {'NaN','nan'}
                    v = NaN;
                otherwise
                    v = str2double(s);
                    if isnan(v)
                        v = 0;
                    end
            end
        end

        %% ---- non-average analyses ------------------------------------------

        function payload = analysisViaCpp(solverName, model, options, analysis, flags)
            % PAYLOAD = ANALYSISVIACPP(SOLVERNAME, MODEL, OPTIONS, ANALYSIS, FLAGS)
            % Run one non-average line-cli analysis and return ITS payload.
            %
            % Every line-cli answer is keyed by its own -a, so reading the key
            % back is also the check that the process answered the question that
            % was put: a payload arriving under another key is a wrong answer,
            % and taking whatever object came would turn that into wrong numbers
            % instead of an error.
            if nargin < 5
                flags = {};
            end
            results = CPPLINE.solveViaCpp(solverName, model, options, analysis, flags);
            key = char(analysis);
            if ~isstruct(results) || ~isfield(results, key)
                keys = {};
                if isstruct(results)
                    keys = fieldnames(results)';
                end
                line_error(mfilename, sprintf('line-cli did not answer -a %s (keys: %s).', ...
                    key, strjoin(keys, ', ')));
            end
            payload = results.(key);
            if ~isstruct(payload)
                line_error(mfilename, sprintf('line-cli''s -a %s payload is not an object.', key));
            end
        end

        function p = probAggr(solverName, model, options, flags)
            % P = PROBAGGR(SOLVERNAME, MODEL, OPTIONS, FLAGS)
            % The -a prob payload as a struct with the fields the four
            % probability getters read: ProbSys, ProbSysAggr and the per-station
            % Prob / ProbAggr vectors, in station order.
            %
            % ALL FOUR COME OFF ONE SOLVE, which is why they share one helper:
            % the CLI computes the stationary law once and reports every view of
            % it, so a getter that ran its own process would pay for the chain
            % again to read another column of the same answer.
            %
            % `-s mva` emits no `ProbSys` / `Prob`: the binomial fit of Schmidt
            % (1997) is an AGGREGATE law with no detailed counterpart, so those
            % fields come back empty and the detailed getters refuse rather than
            % quoting the aggregate. `-s nc` does emit both, from solver_nc_marg
            % and solver_nc_joint, which carry the class-within-chain split the
            % aggregate pair sums out.
            % FLAGS is `getProb`'s second argument on the wire: --node and
            % --state name one node's encoded row, which the arms that accept it
            % substitute for that node's declared state before evaluating.
            if nargin < 4
                flags = {};
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'prob', flags);
            p = struct('ProbSys', [], 'ProbSysAggr', [], 'Prob', [], 'ProbAggr', []);
            if isfield(payload, 'ProbSys')
                p.ProbSys = CPPLINE.jsonNumericScalar(payload.ProbSys);
            end
            if isfield(payload, 'ProbSysAggr')
                p.ProbSysAggr = CPPLINE.jsonNumericScalar(payload.ProbSysAggr);
            end
            if isfield(payload, 'Prob')
                p.Prob = CPPLINE.jsonNumericList(payload.Prob);
            end
            if isfield(payload, 'ProbAggr')
                p.ProbAggr = CPPLINE.jsonNumericList(payload.ProbAggr);
            end
        end

        function v = probEntry(p, field, ist, solverName, method)
            % V = PROBENTRY(P, FIELD, IST, SOLVERNAME, METHOD)
            % One station's entry of a -a prob vector, or a named refusal when
            % the solver's arm does not report that view at all.
            v = [];
            if isfield(p, field)
                v = p.(field);
            end
            if isempty(v)
                CPPLINE.cppUnsupported(solverName, method, sprintf(['line-cli''s -a prob ' ...
                    'answer for this solver carries no %s column: the aggregate law it ' ...
                    'reports has no detailed counterpart'], field));
            end
            if ist < 1 || ist > numel(v)
                line_error(mfilename, sprintf('station index %d is outside 1..%d.', ist, numel(v)));
            end
            v = v(ist);
        end

        function [Pi_t, SS] = tranProb(solverName, model, options, nodeIndex, aggregate)
            % [PI_T, SS] = TRANPROB(SOLVERNAME, MODEL, OPTIONS, NODEINDEX, AGGREGATE)
            % The four getTranProb* answers from -a tranprob: PI_T is [t, pi(t)]
            % and SS is the labelling of the columns of pi(t).
            %
            % PI_T IS THE SAME OBJECT IN ALL FOUR. Every one of the reference
            % getters takes the transient analyzer's DETAILED law and differs
            % only in which labelling it returns beside it, so the aggregate
            % ones read `labelsAggr` and NOT `pitAggr` -- that is the law over
            % aggregate states, a different quantity that would sum to one just
            % as convincingly under the wrong name.
            %
            % SCOPE IS THE --node FLAG: without it the labels span the system,
            % with it they are that node's block, which is the slice the node
            % getters cut out of SS themselves.
            tspan = [];
            if isstruct(options) && isfield(options, 'timespan')
                tspan = options.timespan;
            end
            if numel(tspan) ~= 2 || ~isfinite(tspan(2))
                line_error(mfilename, sprintf(['getTranProb* in %s requires a finite ' ...
                    'timespan T, e.g. %s(model, ''timespan'', [0, T]).'], ...
                    char(solverName), char(solverName)));
            end
            flags = {'--tspan', sprintf('%.17g:%.17g', tspan(1), tspan(2))};
            if ~isempty(nodeIndex)
                flags = [flags, {'--node', sprintf('%d', round(nodeIndex))}];
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'tranprob', flags);
            t = CPPLINE.jsonNumericList(payload.t)';
            Pi_t = [t, CPPLINE.rowsToMatrix(payload.pit)];
            if aggregate
                SS = CPPLINE.rowsToMatrix(payload.labelsAggr);
            else
                SS = CPPLINE.rowsToMatrix(payload.labels);
            end
        end

        function S = samplePath(solverName, model, options, numEvents, nodeIndex)
            % S = SAMPLEPATH(SOLVERNAME, MODEL, OPTIONS, NUMEVENTS, NODEINDEX)
            % The -a sample payload, decoded into the pieces the four sample
            % getters read: the epoch of each step, the system state in both
            % views, the node's own two views, and the synchronization that
            % fired.
            %
            % ONE WALK ANSWERS ALL FOUR, as with probAggr. A trajectory is not a
            % mean: two walks of the same chain are two different answers, so a
            % getter that ran its own process would report a trace that never
            % lines up with the one its sibling reported.
            %
            % THE STATE SPACE COMES BACK WITH IT (`space`, `NodeWidths`), and
            % the visited rows are indices into THAT enumeration. Resolving them
            % against a state space built here instead would pair two orderings
            % that have no reason to agree, and the mismatch would read as a
            % trajectory rather than as an error.
            if nargin < 5
                nodeIndex = [];
            end
            flags = {'--events', sprintf('%d', round(numEvents))};
            % The seed reaches this analysis ONLY here: optionFlags forwards it
            % for the simulator tokens, and line-cli refuses it on every CTMC
            % analysis but this one, so it cannot be forwarded by token.
            if isstruct(options) && isfield(options, 'seed') && ...
                    ~isempty(options.seed) && options.seed > 0
                flags = [flags, {'--seed', sprintf('%d', round(options.seed))}];
            end
            if ~isempty(nodeIndex)
                flags = [flags, {'--node', sprintf('%d', round(nodeIndex))}];
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'sample', flags);
            S = struct('t', [], 'sysState', [], 'nodeWidths', [], ...
                'sysAggr', [], 'nodeState', [], 'nodeAggr', [], 'event', []);
            if isfield(payload, 't')
                S.t = CPPLINE.jsonNumericList(payload.t)';
            end
            if isfield(payload, 'NodeWidths')
                S.nodeWidths = CPPLINE.jsonNumericList(payload.NodeWidths);
            end
            % `state` IS THE TRAJECTORY IN BOTH SPELLINGS. The CTMC sampler
            % walks an enumerated chain and names the rows it visited, so its
            % `state` is a vector of indices into the `space` that travels with
            % it; the SSA sampler has no enumeration to index and sends the
            % visited encodings themselves. Either way what this returns is the
            % system state at each epoch.
            if isfield(payload, 'space')
                space = CPPLINE.rowsToMatrix(payload.space);
                idx = CPPLINE.jsonNumericList(payload.state) + 1;  % indexBase 0
                S.sysState = space(idx, :);
            elseif isfield(payload, 'state')
                S.sysState = CPPLINE.rowsToMatrix(payload.state);
            end
            if isfield(payload, 'sysAggr')
                S.sysAggr = CPPLINE.rowsToMatrix(payload.sysAggr);
            end
            if isfield(payload, 'nodeState')
                S.nodeState = CPPLINE.rowsToMatrix(payload.nodeState);
            end
            if isfield(payload, 'nodeAggr')
                S.nodeAggr = CPPLINE.rowsToMatrix(payload.nodeAggr);
            end
            if isfield(payload, 'event')
                S.event = CPPLINE.sampleEventIndices(payload.event);
            end
        end

        function S = nodeSamplePath(solverName, model, options, node, numEvents, aggregate)
            % S = NODESAMPLEPATH(SOLVERNAME, MODEL, OPTIONS, NODE, NUMEVENTS, AGGREGATE)
            % One node's sample path in the contract of @@SolverCTMC/sample and
            % sampleAggr: handle, t, state, event, isaggregate.
            %
            % BOTH VIEWS COME OFF THE WIRE (`nodeState`, `nodeAggr`). The
            % aggregate one is the port's own marginal_of rather than a
            % State.toMarginal applied here, so the reported counts are the
            % engine's reading of its own states.
            ind = node;
            if isa(node, 'Node')
                ind = node.index;
            end
            P = CPPLINE.samplePath(solverName, model, options, numEvents, ind);
            sn = model.getStruct();
            S = struct();
            S.handle = node;
            S.t = P.t(:);
            if aggregate
                S.state = P.nodeAggr;
            else
                S.state = P.nodeState;
            end
            S.event = CPPLINE.syncEvents(sn, P.event, S.t);
            S.isaggregate = aggregate;
        end

        function S = sysSamplePath(solverName, model, options, numEvents, aggregate)
            % S = SYSSAMPLEPATH(SOLVERNAME, MODEL, OPTIONS, NUMEVENTS, AGGREGATE)
            % The system sample path in the contract of @@SolverCTMC/sampleSys
            % and sampleSysAggr: state is a CELL PER STATEFUL NODE.
            %
            % The blocks are cut out of the wire's own `space` with its own
            % `NodeWidths` -- the local-space widths, which is the padding rule
            % sampleSys uses and sample.m does not.
            %
            % THE AGGREGATE IS TAKEN HERE, with State.toMarginal, and not read
            % off the payload's `sysAggr`: that matrix is (nstations x nclasses)
            % in station order with Sources dropped, while this getter's
            % contract is one entry per STATEFUL node. The two index different
            % things, so quoting one under the other's name would misattribute
            % whole columns. What the engine decided -- which state, at what
            % epoch -- is what P carries; only the aggregation is done here, by
            % the same function the native getter calls.
            P = CPPLINE.samplePath(solverName, model, options, numEvents);
            sn = model.getStruct();
            S = struct();
            S.handle = model.getStatefulNodes';
            S.t = P.t(:);
            nst = cumsum([1, P.nodeWidths]);
            S.state = cell(1, numel(P.nodeWidths));
            for isf = 1:numel(P.nodeWidths)
                block = P.sysState(:, nst(isf):(nst(isf+1)-1));
                if aggregate
                    [~, block] = State.toMarginal(sn, sn.statefulToNode(isf), block);
                end
                S.state{isf} = block;
            end
            S.event = CPPLINE.syncEvents(sn, P.event, S.t);
            S.isaggregate = aggregate;
        end

        function idx = sampleEventIndices(v)
            % IDX = SAMPLEEVENTINDICES(V)
            % The `event` column of -a sample as one-based synchronization
            % indices, with NaN where the walk absorbed.
            %
            % NOT jsonNumericList: that maps a JSON null to 0, and a 0 here
            % becomes index 1 after the base shift -- naming a synchronization
            % that did not fire, at the one step where none did.
            if iscell(v)
                idx = nan(1, numel(v));
                for k = 1:numel(v)
                    if ~isempty(v{k}) && isnumeric(v{k})
                        idx(k) = double(v{k}(1)) + 1;
                    end
                end
            elseif isempty(v)
                idx = [];
            else
                idx = double(v(:)') + 1;
            end
        end

        function E = syncEvents(sn, eventIdx, t)
            % E = SYNCEVENTS(SN, EVENTIDX, T)
            % The sample getters' `event` cell: every active and passive arm of
            % each fired synchronization, stamped with the epoch it fired at.
            %
            % The synchronization DESCRIPTORS are model structure and are read
            % from sn, as the native getters read them; what the C++ decided is
            % WHICH one fired and WHEN, and that is what eventIdx and t carry.
            E = {};
            for e = 1:numel(eventIdx)
                if ~isfinite(eventIdx(e)) || eventIdx(e) < 1 || eventIdx(e) > numel(sn.sync)
                    continue
                end
                s = sn.sync{eventIdx(e)};
                for a = 1:numel(s.active)
                    E{end+1} = s.active{a}; %#ok<AGROW>
                    E{end}.t = t(e);
                end
                for p = 1:numel(s.passive)
                    E{end+1} = s.passive{p}; %#ok<AGROW>
                    E{end}.t = t(e);
                end
            end
        end

        function flags = stateFlags(nodeIndex, state)
            % FLAGS = STATEFLAGS(NODEINDEX, STATE)
            % `getProb(node, state)`'s second argument as argv: the node it
            % belongs to and its ENCODED ROW, comma separated. Empty when no
            % state was passed, which leaves the model's declared one standing.
            flags = {};
            if isempty(state)
                return
            end
            flags = {'--node', sprintf('%d', round(nodeIndex)), '--state', ...
                strjoin(arrayfun(@(v) sprintf('%d', round(v)), state(:)', ...
                'UniformOutput', false), ',')};
        end

        function [Pn, lPn] = probSysMarg(solverName, model, options, nvec, engine)
            % [PN, LPN] = PROBSYSMARG(SOLVERNAME, MODEL, OPTIONS, NVEC, ENGINE)
            % getProbSysMarg: the joint law of the per-station TOTALS, read at
            % NVEC out of the whole lattice -a sysmarg reports.
            %
            % LPN IS log(PN) HERE and not an independently accumulated log: the
            % port's solver_nc_jointmarg returns the probability alone, where
            % the reference also carries the logarithm through. On a population
            % whose probability underflows, this reports -Inf where the native
            % getter still has a finite log.
            sn = model.getStruct();
            if numel(nvec) ~= sn.nstations
                line_error(mfilename, sprintf(['getProbSysMarg takes one job count per ' ...
                    'station: %d given for %d stations.'], numel(nvec), sn.nstations));
            end
            flags = {};
            if nargin >= 5 && ~isempty(engine) && ~strcmpi(char(engine), 'exact')
                flags = {'--perm-engine', lower(char(engine))};
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'sysmarg', flags);
            entries = payload.states;
            if isstruct(entries)
                entries = num2cell(entries);
            end
            want = round(nvec(:)');
            for k = 1:numel(entries)
                if isequal(round(CPPLINE.jsonNumericList(entries{k}.n)), want)
                    Pn = CPPLINE.jsonNumericScalar(entries{k}.P);
                    lPn = log(Pn);
                    return
                end
            end
            line_error(mfilename, sprintf(['line-cli''s -a sysmarg lattice carries no state ' ...
                '[%s]; it sweeps every per-station total summing to the closed population.'], ...
                strjoin(arrayfun(@(v) sprintf('%d', v), want, 'UniformOutput', false), ' ')));
        end

        function [joint, marginal] = mamProb(solverName, model, options, nodeIndex)
            % [JOINT, MARGINAL] = MAMPROB(SOLVERNAME, MODEL, OPTIONS, NODEINDEX)
            % The queue-length law of ONE node under -s mam -a prob: JOINT is
            % the (level x phase) table getProb returns and MARGINAL{r} is the
            % per-class vector getProbMarg returns.
            %
            % BOTH COME OFF ONE SOLVE, because they are two views of the same
            % QBD law and the arm reports them together. The multi-queue refusal
            % is the REFERENCE'S, raised there rather than here, so the two
            % backends decline the same models for the same stated reason.
            flags = {'--node', sprintf('%d', round(nodeIndex))};
            % An open model's level axis is unbounded and the table has to stop
            % somewhere; the cutoff is the caller's, as it is natively.
            if isstruct(options) && isfield(options, 'cutoff') && ~isempty(options.cutoff)
                c = max(options.cutoff(:));
                if isfinite(c) && c > 0
                    flags = [flags, {'--cutoff', CPPLINE.num2arg(c)}];
                end
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'prob', flags);
            joint = CPPLINE.rowsToMatrix(payload.joint);
            entries = payload.marginal;
            if isstruct(entries)
                entries = num2cell(entries);
            end
            marginal = cell(1, numel(entries));
            for k = 1:numel(entries)
                marginal{double(entries{k}.jobclass) + 1} = ...
                    CPPLINE.jsonNumericList(entries{k}.P);
            end
        end

        function lG = normConstAggr(solverName, model, options)
            % LG = NORMCONSTAGGR(SOLVERNAME, MODEL, OPTIONS)
            % getProbNormConstAggr: log G, from -a normconst.
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'normconst');
            if ~isfield(payload, 'logNormConstAggr')
                line_error(mfilename, 'line-cli''s -a normconst answer carries no logNormConstAggr.');
            end
            lG = CPPLINE.jsonNumericScalar(payload.logNormConstAggr);
        end

        function [Pmarg, logPmarg] = probMarg(solverName, model, options, ist, jobclass, state_m)
            % [PMARG, LOGPMARG] = PROBMARG(SOLVERNAME, MODEL, OPTIONS, IST, JOBCLASS, STATE_M)
            % getProbMarg's curve from -a marg, for the station IST.
            %
            % THE TWO SOLVERS ANSWER DIFFERENT QUESTIONS UNDER ONE NAME, which
            % is why the payload key differs and this reads both: SolverMVA's
            % getProbMarg is P(n jobs OF CLASS r), so its curves are keyed by
            % (station, class) and it takes --class and --marg-states; SolverNC's
            % is P(n jobs IN TOTAL at the station), so its curves are keyed by
            % station alone and neither flag applies. Reporting one under the
            % other's name would swap a per-class law for a total-occupancy one.
            sn = model.getStruct();
            if ist < 1 || ist > sn.nstations
                line_error(mfilename, sprintf('station index %d is outside 1..%d.', ist, sn.nstations));
            end
            flags = {'--node', sprintf('%d', sn.stationToNode(ist))};
            if nargin >= 5 && ~isempty(jobclass)
                flags = [flags, {'--class', sprintf('%d', round(jobclass))}];
            end
            if nargin >= 6 && ~isempty(state_m)
                flags = [flags, {'--marg-states', strjoin(arrayfun(@(v) sprintf('%d', ...
                    round(v)), state_m(:)', 'UniformOutput', false), ',')}];
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'marg', flags);
            entries = {};
            if isfield(payload, 'marginal')
                entries = payload.marginal;
            elseif isfield(payload, 'curves')
                entries = payload.curves;
            end
            if isstruct(entries)
                entries = num2cell(entries);
            end
            Pmarg = [];
            logPmarg = [];
            for k = 1:numel(entries)
                e = entries{k};
                % `station` is 0-based on the wire; --node already narrowed the
                % answer, so the first entry is the one asked for, but the index
                % is checked rather than assumed.
                if double(e.station) + 1 ~= ist
                    continue
                end
                Pmarg = CPPLINE.jsonNumericList(e.P);
                logPmarg = CPPLINE.jsonNumericList(e.logP);
                return
            end
            line_error(mfilename, sprintf(['line-cli''s -a marg answer carries no curve for ' ...
                'station %d.'], ist));
        end

        % REFUSESTATEDEPENDENT IS GONE, and its absence is the point. It stood
        % for "this getter has no line-cli analysis to delegate to", and its 23
        % call sites across SolverCTMC, SolverSSA, SolverMVA, SolverNC and
        % SolverMAM were checked one by one against what the CLI dispatches
        % today: every one of them had an arm waiting (-a sample, -a tranprob,
        % -a prob, -a marg, -a sysmarg, -a normconst). A blanket refusal helper
        % is what let those claims decay unnoticed, so a getter that genuinely
        % cannot be served now calls cppUnsupported with ITS OWN reason, which
        % names the thing that is missing and can therefore be rechecked.

        function assertSingleState(solverName, method, model)
            % ASSERTSINGLESTATE(SOLVERNAME, METHOD, MODEL)
            % Refuse a state-dependent analysis whose initial state the wire
            % cannot carry.
            %
            % THE SINGLE STATE NOW TRAVELS. linemodel_save emits the
            % (stateSpace, statePrior) pair for every stateful node carrying a
            % state, INCLUDING the one-row space with prior [1] that setState
            % and initFromMarginal leave, and the C++ network_reader stores it;
            % analyzer_detail::default_init_state takes that row in preference
            % to the default marking, and api::sn_declared_marginal decodes it
            % for the -a prob arms. So -a prob, the transient arms and a sampled
            % trajectory all start where the caller put the jobs, and this used
            % to refuse them all -- unconditionally, because MATLAB has no
            % marker telling a setState apart from an initDefault. That marker
            % is no longer needed: the two are now sent alike, and where the
            % state IS the default the C++ reconstructs the same row.
            %
            % WHAT IS STILL REFUSED IS A GENUINE DISTRIBUTION over several rows,
            % AND ONLY ON THE ARMS THAT SEED ONE STATE. default_init_state and
            % sn_declared_marginal take row 0 only when the space has exactly
            % ONE row and rebuild the default marking otherwise, so a prior over
            % k > 1 states would be answered as one state under the name of a
            % question about a mixture.
            %
            % IT DOES NOT APPLY TO THE TRANSIENT ARMS.
            % solver_ctmc_transient_analyzer seeds init_state_distribution, the
            % product of the declared per-node priors over the enumerated space,
            % and integrates ONCE from that mixture, which IS the reference's
            % weighted sum over the support. -a tran, -a tranprob and
            % -a tranreward therefore take a mixture and are not gated here.
            % See _kb/07-cross-language-parity.md.
            if nargin < 3 || isempty(model)
                return
            end
            named = {};
            nodes = model.getNodes();
            for ni = 1:numel(nodes)
                nd = nodes{ni};
                if ~isa(nd, 'StatefulNode'), continue; end
                if numel(nd.statePrior) > 1 || size(nd.space, 1) > 1
                    named{end+1} = nd.getName(); %#ok<AGROW>
                end
            end
            if isempty(named)
                return
            end
            CPPLINE.cppUnsupported(solverName, method, sprintf(['the initial state of %s is a ' ...
                'distribution over several states: the C++ reader takes a one-row state space ' ...
                'as the initial state and rebuilds the default marking for any other, so it ' ...
                'would answer for one state under a question about a mixture. The reference ' ...
                'weights the analysis over the prior''s support'], strjoin(named, ', ')));
        end

        function RD = cdfRespT(solverName, model, options)
            % RD = CDFRESPT(SOLVERNAME, MODEL, OPTIONS)
            % The per-station, per-class response-time CDF from -a cdf, in
            % @SolverCTMC/getCdfRespT's own contract: an (nstations x nclasses)
            % cell whose entries are [F(t), t] two-column matrices.
            %
            % This is the EXACT tagged-chain law, the same quantity the MATLAB
            % getter builds, so the two are comparable value for value. A
            % (station, class) pair the tagged chain never visits carries no
            % curve on the wire and stays [] here, which is how the MATLAB getter
            % also leaves it.
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'cdf');
            sn = model.getStruct();
            RD = cell(sn.nstations, sn.nclasses);
            if ~isfield(payload, 'respt')
                return
            end
            entries = payload.respt;
            if isstruct(entries)
                entries = num2cell(entries);
            end
            for k = 1:numel(entries)
                e = entries{k};
                % indexBase is 0 on the wire against the tables' 1-based columns
                ist = double(e.station) + 1;
                r = double(e.jobclass) + 1;
                t = double(e.t(:));
                F = double(e.F(:));
                if numel(t) ~= numel(F)
                    line_error(mfilename, sprintf(['line-cli sent %d abscissae and %d CDF ' ...
                        'values for station %d class %d; a curve must have one of each.'], ...
                        numel(t), numel(F), ist, r));
                end
                RD{ist, r} = [F, t];
            end
        end

        function spec = passageSetSpec(S)
            % SPEC = PASSAGESETSPEC(S)
            % A state set as the --passage-from/--passage-into wire syntax:
            % rows joined by ';', entries by ','. STATE ROWS travel rather than
            % row indices, because the two enumerations need not order (or even
            % purge) states identically -- line-cli resolves rows by content
            % against the space its own engine enumerated.
            rows = cell(1, size(S,1));
            for i = 1:size(S,1)
                rows{i} = strjoin(arrayfun(@(v) sprintf('%.17g', v), S(i,:), ...
                    'UniformOutput', false), ',');
            end
            spec = strjoin(rows, ';');
        end

        function [RD, out] = cdfFirstPassT(solverName, model, options, Arows, Brows)
            % [RD, OUT] = CDFFIRSTPASST(SOLVERNAME, MODEL, OPTIONS, AROWS, BROWS)
            % The state-set first passage law from -a firstpasst, in
            % @SolverCTMC/getCdfFirstPassT's own contract: RD an [n x 2] matrix
            % of [F(t), t], OUT carrying tset and the density. AROWS/BROWS are
            % STATE ROWS (already resolved by the caller against its own space);
            % an empty AROWS selects the conditional stationary law.
            flags = {};
            if ~isempty(Arows)
                flags{end+1} = '--passage-from'; %#ok<AGROW>
                flags{end+1} = CPPLINE.passageSetSpec(Arows); %#ok<AGROW>
            end
            flags{end+1} = '--passage-into';
            flags{end+1} = CPPLINE.passageSetSpec(Brows);
            config = options.config;
            if isfield(config,'passage_method') && ~isempty(config.passage_method)
                flags{end+1} = '--passage-method'; %#ok<AGROW>
                flags{end+1} = char(config.passage_method); %#ok<AGROW>
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'firstpasst', flags);
            t = double(payload.t(:));
            F = double(payload.F(:));
            if numel(t) ~= numel(F)
                line_error(mfilename, sprintf(['line-cli sent %d abscissae and %d CDF values ' ...
                    'for the first passage; a curve must have one of each.'], numel(t), numel(F)));
            end
            RD = [F, t];
            out = struct();
            out.tset = t;
            out.density = double(payload.f(:));
        end

        function [m, mall] = firstPassTMoments(solverName, model, options, Arows, Brows, nmax)
            % [M, MALL] = FIRSTPASSTMOMENTS(SOLVERNAME, MODEL, OPTIONS, AROWS, BROWS, NMAX)
            % The state-set first passage MOMENTS from -a firstpasstmom, in
            % @SolverCTMC/getFirstPassTMoments's own contract: M the (1 x NMAX)
            % moment vector, MALL the (nstates x NMAX) per-source matrix.
            % AROWS/BROWS are STATE ROWS, as in CDFFIRSTPASST.
            %
            % No --passage-method here, and that is not an omission: the moments
            % come from one linear solve per order and invert no transform, so
            % there is no inversion for the flag to select. line-cli refuses it
            % on this arm for the same reason.
            flags = {};
            if ~isempty(Arows)
                flags{end+1} = '--passage-from'; %#ok<AGROW>
                flags{end+1} = CPPLINE.passageSetSpec(Arows); %#ok<AGROW>
            end
            flags{end+1} = '--passage-into';
            flags{end+1} = CPPLINE.passageSetSpec(Brows);
            flags{end+1} = '--passage-orders';
            flags{end+1} = num2str(round(nmax));
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'firstpasstmom', flags);
            m = double(payload.m(:))';
            mall = [];
            if isfield(payload, 'mall') && ~isempty(payload.mall)
                rows = payload.mall;
                if iscell(rows)
                    mall = zeros(numel(rows), numel(m));
                    for i = 1:numel(rows)
                        mall(i,:) = double(rows{i}(:))';
                    end
                else
                    mall = double(rows);
                end
            end
        end

        function [J, rhs, vars, equilibria] = jacobian(solverName, model, options, wantEquilibria)
            % [J, RHS, VARS, EQUILIBRIA] = JACOBIAN(SOLVERNAME, MODEL, OPTIONS, WANTEQUILIBRIA)
            % The SYMBOLIC drift Jacobian from -a jacobian, in
            % @SolverFLD/getJacobian's own contract.
            %
            % This is symbolic on both sides, which is why it can be bridged at
            % all: line-cli's arm builds the drift with fluid_symodes and
            % differentiates it through the SAME symbolic backend the MATLAB
            % getter reaches through SAGE, so the entries are expressions and
            % not a numeric evaluation of them.
            %
            % EQUILIBRIA are only computed when asked for -- solving the system
            % is the expensive part -- so WANTEQUILIBRIA carries the CALLER's
            % nargout, not this function's: called as [J,rhs,vars,eq] = ...,
            % nargout here is always 4 and would ask for them every time.
            % hasEquilibria on the wire separates "asked and there are none"
            % from "never asked", so an empty list is never silently reported as
            % "this system has none".
            flags = {};
            if nargin >= 4 && wantEquilibria
                flags = {'--equilibria'};
            end
            options = CPPLINE.pnormMethod(options);
            p = CPPLINE.analysisViaCpp(solverName, model, options, 'jacobian', flags);
            J = {};
            if isfield(p, 'jacobian')
                J = SAGE.asCellMatrix(p.jacobian);
            end
            rhs = {};
            if isfield(p, 'rhs')
                rhs = CPPLINE.jsonStringList(p.rhs);
            end
            vars = {};
            if isfield(p, 'vars')
                vars = CPPLINE.jsonStringList(p.vars);
            end
            equilibria = {};
            if isfield(p, 'equilibria')
                equilibria = p.equilibria;
            end
        end

        function options = pnormMethod(options)
            % OPTIONS = PNORMMETHOD(OPTIONS)
            % Restate MATLAB's p-norm smoothing request in line-cli's own terms.
            %
            % THE TWO SIDES SELECT THE SMOOTHING DIFFERENTLY, and only the
            % symbolic getters notice, because a min-scaled drift has no
            % Jacobian and so must be smoothed before it can be differentiated
            % at all. MATLAB turns the smoothing on whenever options.config.pstar
            % is set, leaving the METHOD name alone (solver_fluid_symodes,
            % use_pnorm); line-cli turns it on only under --method pnorm and has
            % no flag carrying the exponent, so it always smooths at its own
            % p = 20. Forwarding 'matrix' would therefore have line-cli refuse a
            % model MATLAB answers for, under a message about a kink the caller
            % had already smoothed away.
            %
            % A DIFFERENT EXPONENT IS REFUSED RATHER THAN ROUNDED TO 20: p sets
            % how sharply the smoothed min approaches the kink, so answering at
            % another p is answering about another drift.
            if isempty(options) || ~isstruct(options)
                return
            end
            pstar = [];
            if isfield(options, 'pstar') && ~isempty(options.pstar)
                pstar = options.pstar;
            elseif isfield(options, 'config') && isstruct(options.config) && ...
                    isfield(options.config, 'pstar') && ~isempty(options.config.pstar)
                pstar = options.config.pstar;
            end
            if isempty(pstar)
                return
            end
            p = unique(double(pstar(:)));
            if ~isscalar(p) || p ~= CPPLINE.CLI_PSTAR
                line_error(mfilename, sprintf(['lang=''cpp'' cannot serve the symbolic ' ...
                    'getters at this smoothing: line-cli selects the p-norm smoothing with ' ...
                    '--method pnorm and carries no flag for its exponent, so it can only ' ...
                    'answer at p = %g, while options.config.pstar asks for %s. Use p = %g, ' ...
                    'or solve with lang=''matlab''.'], CPPLINE.CLI_PSTAR, ...
                    mat2str(p(:)'), CPPLINE.CLI_PSTAR));
            end
            options.method = 'pnorm';
        end

        function [t, QVart, Sigmat] = tranAvgVar(solverName, model, options)
            % [T, QVART, SIGMAT] = TRANAVGVAR(SOLVERNAME, MODEL, OPTIONS)
            % The transient second moments from -a tranvar, in
            % @SolverFLD/getTranAvgVar's contract: T the abscissae, QVART an
            % (nstations x nclasses) cell of variance series, and SIGMAT the
            % (dim x dim x numel(t)) covariance trajectory.
            %
            % THE HORIZON IS ALREADY RESOLVED by the caller, which applies the
            % same slowest-rate rule the native path uses for an unbounded
            % timespan; this only transmits it. line-cli enforces the kp
            % precondition on its own side, so both agree that the other fluid
            % methods integrate the mean alone.
            tspan = options.timespan;
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'tranvar', ...
                {'--tspan', sprintf('%.17g:%.17g', tspan(1), tspan(2))});
            sn = model.getStruct();
            t = [];
            if isfield(payload, 't')
                t = CPPLINE.jsonNumericList(payload.t)';
            end
            QVart = cell(sn.nstations, sn.nclasses);
            if isfield(payload, 'curves')
                entries = payload.curves;
                if isstruct(entries)
                    entries = num2cell(entries);
                end
                for k = 1:numel(entries)
                    e = entries{k};
                    QVart{double(e.station) + 1, double(e.jobclass) + 1} = ...
                        CPPLINE.jsonNumericList(e.QVar)';
                end
            end
            % Sigma arrives as one matrix PER TIME POINT; the getter's contract
            % is (dim, dim, numel(t)), so TIME HAS TO BECOME THE THIRD INDEX.
            % jsondecode collapses the rectangular stack into a numeric array
            % whose FIRST index is time, and num2cell(S,[1 2]) then slices it
            % along the wrong pair: it cuts the (nt x dim) planes and returned a
            % (nt x dim x dim) array, so a caller reading Sigmat(:,:,k) got a
            % time series where a covariance matrix was promised.
            Sigmat = [];
            if isfield(payload, 'Sigma')
                S = payload.Sigma;
                if iscell(S)
                    for n = 1:numel(S)
                        B = CPPLINE.rowsToMatrix(S{n});
                        if n == 1
                            Sigmat = zeros(size(B, 1), size(B, 2), numel(S));
                        end
                        Sigmat(:, :, n) = B; %#ok<AGROW>
                    end
                else
                    A = double(S);
                    nt = numel(t);
                    if nt > 0 && size(A, 1) == nt && mod(numel(A), nt) == 0
                        % dim is recovered from the element count because
                        % jsondecode drops a trailing singleton dimension, which
                        % is what a one-state system sends.
                        dim = round(sqrt(numel(A) / nt));
                        Sigmat = permute(reshape(A, nt, dim, dim), [2 3 1]);
                    else
                        % A single time point arrives already squeezed to the
                        % (dim x dim) block, which is (dim x dim x 1).
                        Sigmat = A;
                    end
                end
            end
        end

        function aoiResults = aoiResults(solverName, model, options)
            % AOIRESULTS = AOIRESULTS(SOLVERNAME, MODEL, OPTIONS)
            % The -a aoi payload in solver_mfq_aoi's OWN struct shape, the one
            % @SolverFLD/getAvgAoI and getCdfAoI already read out of
            % self.result.solverSpecific.aoiResults.
            %
            % Returning THAT struct rather than the finished answers is what
            % keeps the port whole: both getters then run their existing bodies,
            % so the summary table and the survival-function algebra have one
            % implementation instead of a MATLAB one and a bridge one that can
            % round differently. The (g,A,h) triples ride along for exactly this
            % reason -- getCdfAoI takes a caller's own t_values, and a fixed grid
            % could not answer that.
            p = CPPLINE.analysisViaCpp(solverName, model, options, 'aoi');
            aoiResults = struct();
            aoiResults.AoI_mean = CPPLINE.jsonScalar(p, 'AoIMean');
            aoiResults.AoI_var = CPPLINE.jsonScalar(p, 'AoIVar');
            aoiResults.PAoI_mean = CPPLINE.jsonScalar(p, 'PAoIMean');
            aoiResults.PAoI_var = CPPLINE.jsonScalar(p, 'PAoIVar');
            aoiResults.systemType = p.systemType;
            aoiResults.preemption = CPPLINE.jsonScalar(p, 'preemption');
            % g is a ROW and h a COLUMN: getCdfAoI evaluates
            % g * expm(A t) / A * h, so the orientations are load-bearing and
            % not a style choice. jsonNumericList answers a row either way.
            aoiResults.AoI_g = CPPLINE.jsonNumericList(p.AoI_g);
            aoiResults.AoI_h = CPPLINE.jsonNumericList(p.AoI_h)';
            aoiResults.AoI_A = CPPLINE.rowsToMatrix(p.AoI_A);
            aoiResults.PAoI_g = CPPLINE.jsonNumericList(p.PAoI_g);
            aoiResults.PAoI_h = CPPLINE.jsonNumericList(p.PAoI_h)';
            aoiResults.PAoI_A = CPPLINE.rowsToMatrix(p.PAoI_A);
        end

        function M = rowsToMatrix(v)
            % M = ROWSTOMATRIX(V)
            % A matrix sent as an array of row arrays. jsondecode already
            % collapses a rectangular one into a numeric matrix, so this only
            % has to reassemble the cell form and keep the empty case 0x0.
            if isempty(v)
                M = [];
            elseif iscell(v)
                rows = cellfun(@(r) CPPLINE.jsonNumericList(r), v(:)', 'UniformOutput', false);
                M = vertcat(rows{:});
            else
                M = double(v);
            end
        end

        function [QNt, UNt, TNt] = tranAvg(solverName, model, options)
            % [QNT, UNT, TNT] = TRANAVG(SOLVERNAME, MODEL, OPTIONS)
            % The transient means from -a tran, in the contract of
            % self.result.Tran.Avg: (nstations x nclasses) cells whose entries
            % are [value, t] two-column matrices.
            %
            % THAT IS THE RESULT TABLE, NOT THE GETTER'S ANSWER. getTranAvg
            % returns a metricVal STRUCT per cell (handle, t, metric,
            % isaggregate) and NaN where a handle is disabled, and it builds
            % those in @NetworkSolver/getTranAvg from exactly these tables. The
            % callers therefore store what comes back and delegate there, so the
            % wrapping has one implementation rather than a native one and a
            % bridge one that can disagree on the disabled case.
            %
            % THE HORIZON IS THE CALLER'S, not this bridge's. line-cli refuses
            % -a tran without --tspan rather than inventing one, so the timespan
            % is read from options here and a model asked for a transient
            % without one is refused on this side, where the option can be named.
            tspan = [];
            if isfield(options, 'timespan')
                tspan = options.timespan;
            end
            if numel(tspan) ~= 2 || ~all(isfinite(tspan)) || ~(tspan(1) < tspan(2))
                line_error(mfilename, ['a transient analysis needs a finite horizon: set ' ...
                    'options.timespan = [t0 t1] with t0 < t1 (for example ' ...
                    'SolverFluid(model, ''timespan'', [0 50])).']);
            end
            flags = {'--tspan', sprintf('%.17g:%.17g', tspan(1), tspan(2))};
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'tran', flags);
            sn = model.getStruct();
            QNt = cell(sn.nstations, sn.nclasses);
            UNt = cell(sn.nstations, sn.nclasses);
            TNt = cell(sn.nstations, sn.nclasses);
            if ~isfield(payload, 'curves')
                return
            end
            entries = payload.curves;
            if isstruct(entries)
                entries = num2cell(entries);
            end
            for k = 1:numel(entries)
                e = entries{k};
                ist = double(e.station) + 1;  % indexBase is 0 on the wire
                r = double(e.jobclass) + 1;
                t = double(e.t(:));
                QNt{ist, r} = [double(e.QLen(:)), t];
                UNt{ist, r} = [double(e.Util(:)), t];
                TNt{ist, r} = [double(e.Tput(:)), t];
            end
        end

        function RD = mamCdfRespT(solverName, model, options)
            % RD = MAMCDFRESPT(SOLVERNAME, MODEL, OPTIONS)
            % The MAM response-time CDF from -a cdf, in solver_mam_passage_time's
            % own contract: an (nstations x nclasses) cell whose QUEUE row holds
            % [F(t), t] and whose Source row stays [].
            %
            % The wire entries carry a class and NO station, and that is not an
            % omission: solver_mam_passage_time covers exactly two stations, a
            % Source and one queue, so the station is determined by the model
            % rather than by the class. This resolves it the same way, and errors
            % if the model is not of that shape rather than guessing a row.
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'cdf');
            sn = model.getStruct();
            RD = cell(sn.nstations, sn.nclasses);
            if ~isfield(payload, 'respt')
                return
            end
            % solver_mam_passage_time's OWN rule for idx_q, not a re-derivation:
            % it takes the queue as the station scheduling FCFS, HOL, FCFSPRPRIO
            % or PS, and the arrival as the EXT one.
            queueRow = find(ismember(sn.sched, [SchedStrategy.FCFS, SchedStrategy.HOL, ...
                SchedStrategy.FCFSPRPRIO, SchedStrategy.PS]));
            if numel(queueRow) ~= 1
                line_error(mfilename, sprintf(['the C++ MAM response-time law is defined for a ' ...
                    'Source and ONE queue, and this model has %d stations scheduling FCFS, HOL, ' ...
                    'FCFSPRPRIO or PS, so the curve line-cli sent could not be placed.'], ...
                    numel(queueRow)));
            end
            entries = payload.respt;
            if isstruct(entries)
                entries = num2cell(entries);
            end
            for k = 1:numel(entries)
                e = entries{k};
                r = double(e.jobclass) + 1;  % indexBase is 0 on the wire
                t = double(e.t(:));
                F = double(e.F(:));
                if numel(t) ~= numel(F)
                    line_error(mfilename, sprintf(['line-cli sent %d abscissae and %d CDF ' ...
                        'values for class %d; a curve must have one of each.'], ...
                        numel(t), numel(F), r));
                end
                RD{queueRow, r} = [F, t];
            end
        end

        function RD = cdfSysRespT(solverName, model, options)
            % RD = CDFSYSRESPT(SOLVERNAME, MODEL, OPTIONS)
            % The PER-CHAIN system response-time CDF, in @SolverCTMC's
            % getCdfSysRespT contract: a (1 x nchains) cell of [F(t), t]. The
            % ROW shape is that getter's, not a convention: it builds
            % cell(1, nchains) and callers index it linearly.
            %
            % It rides in the SAME -a cdf payload as the per-station curves,
            % under 'sysrespt', because line-cli's solve_ctmc_cdf computes both
            % from one tagged-chain solve. Asking for it separately would solve
            % the chain twice and report the second one.
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'cdf');
            sn = model.getStruct();
            RD = cell(1, sn.nchains);
            if ~isfield(payload, 'sysrespt')
                return
            end
            entries = payload.sysrespt;
            if isstruct(entries)
                entries = num2cell(entries);
            end
            for k = 1:numel(entries)
                e = entries{k};
                c = double(e.chain) + 1;  % indexBase is 0 on the wire
                t = double(e.t(:));
                F = double(e.F(:));
                if numel(t) ~= numel(F)
                    line_error(mfilename, sprintf(['line-cli sent %d abscissae and %d CDF ' ...
                        'values for chain %d; a curve must have one of each.'], ...
                        numel(t), numel(F), c));
                end
                RD{c} = [F, t];
            end
        end

        function [R, names] = avgReward(solverName, model, options)
            % [R, NAMES] = AVGREWARD(SOLVERNAME, MODEL, OPTIONS)
            % The steady-state expectation of every declared reward, from
            % -a reward, in @SolverCTMC/getAvgReward's contract.
            %
            % The serializability precondition is shared with getTranReward;
            % see unbridgeableRewards for why a dropped reward must stop the
            % -a reward route rather than shorten its answer, and rewardsLocally
            % for the route it takes instead.
            if ~isempty(CPPLINE.unbridgeableRewards(model))
                [R, names] = CPPLINE.avgRewardLocally(solverName, model, options);
                return
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'reward');
            R = [];
            names = {};
            if isfield(payload, 'E')
                R = CPPLINE.jsonNumericList(payload.E)';
            end
            if isfield(payload, 'Reward')
                names = CPPLINE.jsonStringList(payload.Reward)';
            end
            [ord, ok] = CPPLINE.rewardWireOrder(model, names);
            if ok
                R = R(ord);
                names = names(ord);
            end
        end

        function [Rt, t, names] = tranReward(solverName, model, options, rewardName)
            % [RT, T, NAMES] = TRANREWARD(SOLVERNAME, MODEL, OPTIONS, REWARDNAME)
            % The transient expectation of every declared reward, from
            % -a tranreward, in @SolverCTMC/getTranReward's contract: RT a
            % (nrewards x 1) cell of time series, T the shared abscissae, NAMES
            % the matching (nrewards x 1) cell. A named reward returns that one
            % alone, as the native getter does.
            %
            % The serializability precondition and the wire-order correction are
            % the SAME as getAvgReward's, and are shared rather than restated:
            % a bare function handle still has no serializable form here, and
            % the wire still sorts by name where the model declares in its own
            % order.
            tspan = [];
            if isfield(options, 'timespan')
                tspan = options.timespan;
            end
            if numel(tspan) ~= 2 || ~all(isfinite(tspan)) || ~(tspan(1) < tspan(2))
                line_error(mfilename, ['a transient reward needs a finite horizon: set ' ...
                    'options.timespan = [t0 t1] with t0 < t1.']);
            end
            if ~isempty(CPPLINE.unbridgeableRewards(model))
                [Rt, t, names] = CPPLINE.tranRewardLocally(solverName, model, options, tspan);
                if nargin >= 4 && ~isempty(rewardName)
                    idx = find(strcmp(names, rewardName), 1);
                    if isempty(idx)
                        line_error(mfilename, sprintf(['reward ''%s'' is not declared on this ' ...
                            'model; the declared rewards are %s.'], rewardName, ...
                            strjoin(names', ', ')));
                    end
                    Rt = Rt{idx};
                    names = names{idx};
                end
                return
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'tranreward', ...
                {'--tspan', sprintf('%.17g:%.17g', tspan(1), tspan(2))});
            t = [];
            if isfield(payload, 't')
                t = CPPLINE.jsonNumericList(payload.t)';
            end
            names = {};
            if isfield(payload, 'Reward')
                names = CPPLINE.jsonStringList(payload.Reward)';
            end
            % EACH ENTRY IS A metricVal STRUCT, not the bare series: the native
            % getter returns .t, .metric and .name per reward and its callers
            % plot Rt{r}.metric against Rt{r}.t, so handing back a column vector
            % would break every one of them at the dot rather than in the
            % numbers. The name travels INSIDE the struct so the wire-order
            % correction below permutes the two together by construction.
            Rt = cell(numel(names), 1);
            if isfield(payload, 'E')
                E = payload.E;
                if ~iscell(E)
                    E = num2cell(E, 2);
                end
                for i = 1:min(numel(E), numel(names))
                    metricVal = struct();
                    metricVal.t = t;
                    metricVal.metric = CPPLINE.jsonNumericList(E{i})';
                    metricVal.name = names{i};
                    Rt{i} = metricVal;
                end
            end
            [ord, ok] = CPPLINE.rewardWireOrder(model, names);
            if ok
                Rt = Rt(ord);
                names = names(ord);
            end
            if nargin >= 4 && ~isempty(rewardName)
                idx = find(strcmp(names, rewardName), 1);
                if isempty(idx)
                    line_error(mfilename, sprintf(['reward ''%s'' is not declared on this ' ...
                        'model; line-cli answered for %s.'], rewardName, strjoin(names', ', ')));
                end
                Rt = Rt{idx};
                names = names{idx};
            end
        end

        function unbridgeable = unbridgeableRewards(model)
            % UNBRIDGEABLE = UNBRIDGEABLEREWARDS(MODEL)
            % The declared rewards that no wire format can carry, by name.
            %
            % A REWARD THE WRITER DROPS MUST STOP THE -a reward ROUTE, NOT
            % SHORTEN ITS ANSWER. linemodel_save omits a reward whose descriptor
            % is absent or Custom -- a bare function handle has no serializable
            % form -- and only warns. Taking that route anyway hands line-cli the
            % SURVIVING rewards, whose shorter vector the caller then pairs with
            % its own full declaration list; where every reward is a handle it
            % hands over none and line-cli exits 2. This list is what sends the
            % solve down rewardsLocally instead, and is shared by the
            % steady-state and the transient arm so the two cannot disagree
            % about which model is bridgeable.
            unbridgeable = {};
            if isprop(model,'sn') && isfield(model.sn,'reward')
                for i = 1:numel(model.sn.reward)
                    rw = model.sn.reward{i};
                    d = [];
                    if isfield(rw,'descriptor'), d = rw.descriptor; end
                    if isempty(d) || ~isa(d,'RewardDescriptor') || strcmp(d.kind,'Custom') ...
                            || isempty(d.node)
                        unbridgeable{end+1} = rw.name; %#ok<AGROW>
                    end
                end
            end
        end

        function [R, names] = rewardMatrixOver(model, spaceAggr)
            % [R, NAMES] = REWARDMATRIXOVER(MODEL, SPACEAGGR)
            % The (nrewards x nstates) matrix of the declared rewards evaluated
            % on the rows of SPACEAGGR, and their names in DECLARATION order.
            %
            % The reward map is a function of the AGGREGATE state row and of
            % nothing else, so the same evaluation serves whichever engine
            % enumerated the space. This is solver_ctmc_reward's inner loop,
            % kept identical down to the two calling conventions it accepts, so
            % a handle answers the same number under lang='cpp' as it does
            % natively -- the only thing that differs is which engine built the
            % space and the law it is weighted against.
            sn = model.getStruct(true);
            nodeToStationMap = configureDictionary('int32', 'int32');
            classToIndexMap = configureDictionary('int32', 'int32');
            for ind = 1:sn.nnodes
                if sn.isstation(ind)
                    nodeToStationMap(int32(ind)) = sn.nodeToStation(ind);
                end
            end
            for r = 1:sn.nclasses
                classToIndexMap(int32(r)) = r;
            end
            nRewards = numel(sn.reward);
            nstates = size(spaceAggr, 1);
            R = zeros(nRewards, nstates);
            names = cell(nRewards, 1);
            for r = 1:nRewards
                names{r} = sn.reward{r}.name;
                rewardFn = sn.reward{r}.fn;
                for s = 1:nstates
                    stateRow = spaceAggr(s, :);
                    rewardState = RewardState(stateRow, sn, nodeToStationMap, classToIndexMap);
                    try
                        R(r, s) = rewardFn(rewardState);
                    catch ME
                        try
                            R(r, s) = rewardFn(stateRow, sn);
                        catch
                            rethrow(ME);
                        end
                    end
                end
            end
        end

        function [R, names] = avgRewardLocally(solverName, model, options)
            % [R, NAMES] = AVGREWARDLOCALLY(SOLVERNAME, MODEL, OPTIONS)
            % E[r] for rewards no wire format can carry, against the C++ law.
            %
            % A FUNCTION HANDLE CANNOT CROSS THE WIRE, and there is nothing to
            % serialize it into -- but the reward is a function OF THE
            % STATIONARY LAW, and that law is line-cli's: pi and the aggregated
            % state space arrive through CPPLINE.generator. Applying the
            % caller's own handle to them therefore still answers under
            % lang='cpp': the engine that computed the chain, its stationary
            % distribution and its enumeration is the C++ one, and this method
            % only evaluates a map that no wire format can express. Refusing
            % instead left rewardModel_aggregation, rewardModel_mm1k and
            % rewardModel_multiclass -- whose rewards are handles -- unsolvable.
            %
            % EVERY reward is evaluated here, not just the unbridgeable ones:
            % -a reward answers for the serializable subset only, and splicing
            % two vectors computed over the same law adds nothing but a pairing
            % that can rotate.
            g = CPPLINE.generator(solverName, model, options);
            if isempty(g.pi) || isempty(g.spaceAggr)
                line_error(mfilename, ['line-cli returned no stationary law for this model, ' ...
                    'so the reward(s) defined by a bare function handle cannot be evaluated ' ...
                    'against it.']);
            end
            if size(g.spaceAggr, 1) ~= numel(g.pi)
                line_error(mfilename, sprintf(['line-cli returned a %d-entry stationary law ' ...
                    'over a %d-row aggregated state space; a reward cannot be paired with ' ...
                    'it.'], numel(g.pi), size(g.spaceAggr, 1)));
            end
            [Rmat, names] = CPPLINE.rewardMatrixOver(model, g.spaceAggr);
            R = Rmat * g.pi(:);
        end

        function [Rt, t, names] = tranRewardLocally(solverName, model, options, tspan)
            % [RT, T, NAMES] = TRANREWARDLOCALLY(SOLVERNAME, MODEL, OPTIONS, TSPAN)
            % E[r(X(t))] for rewards no wire format can carry, against the C++
            % transient law.
            %
            % The twin of avgRewardLocally, over pi(t) instead of pi: -a
            % tranprob integrates the forward equation ONCE from the model's own
            % initial distribution and returns the trajectory beside the
            % AGGREGATED labels its columns are indexed by, which is exactly the
            % pair the reward map needs.
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'tranprob', ...
                {'--tspan', sprintf('%.17g:%.17g', tspan(1), tspan(2))});
            t = [];
            if isfield(payload, 't')
                t = CPPLINE.jsonNumericList(payload.t)';
            end
            pit = [];
            if isfield(payload, 'pitAggr')
                pit = CPPLINE.rowsToMatrix(payload.pitAggr);
            end
            spaceAggr = [];
            if isfield(payload, 'labelsAggr')
                spaceAggr = CPPLINE.rowsToMatrix(payload.labelsAggr);
            end
            if isempty(t) || isempty(pit) || isempty(spaceAggr)
                line_error(mfilename, ['line-cli returned no transient law for this model, ' ...
                    'so the reward(s) defined by a bare function handle cannot be evaluated ' ...
                    'against it.']);
            end
            if size(pit, 2) ~= size(spaceAggr, 1)
                line_error(mfilename, sprintf(['line-cli returned a %d-column transient law ' ...
                    'over a %d-row aggregated state space; a reward cannot be paired with ' ...
                    'it.'], size(pit, 2), size(spaceAggr, 1)));
            end
            [Rmat, names] = CPPLINE.rewardMatrixOver(model, spaceAggr);
            E = pit * Rmat';   % (ntimes x nrewards)
            Rt = cell(numel(names), 1);
            for i = 1:numel(names)
                metricVal = struct();
                metricVal.t = t;
                metricVal.metric = E(:, i);
                metricVal.name = names{i};
                Rt{i} = metricVal;
            end
        end

        function [order, ok] = rewardWireOrder(model, names)
            % [ORDER, OK] = REWARDWIREORDER(MODEL, NAMES)
            % The permutation taking line-cli's answer back into DECLARATION
            % order.
            %
            % THE WIRE SORTS THE REWARDS; THE MODEL DOES NOT. linemodel_save
            % emits them in NAME order so the three writers byte-match, so
            % line-cli answers in that order while the native getter answers in
            % DECLARATION order. A caller pairing the values with its own
            % declaration list then reads every reward under the wrong name --
            % on rewardModel_templates a three-way rotation, not a numeric error.
            order = [];
            ok = false;
            if isempty(names) || ~isprop(model,'sn') || ~isfield(model.sn,'reward')
                return
            end
            if numel(model.sn.reward) ~= numel(names)
                return
            end
            order = zeros(numel(names),1);
            for i = 1:numel(model.sn.reward)
                pos = find(strcmp(names, model.sn.reward{i}.name), 1);
                if isempty(pos)
                    order = [];
                    return
                end
                order(i) = pos;
            end
            ok = true;
        end

        function g = generator(solverName, model, options)
            % G = GENERATOR(SOLVERNAME, MODEL, OPTIONS)
            % The CTMC generator, state spaces and stationary law, as a struct
            % with fields infGen, space, spaceAggr, pi and eventFilt.
            %
            % TWO INVOCATIONS, ONE PER GETTER: -a gen is getGenerator and
            % -a states is getStateSpace, and line-cli keeps them apart because
            % they are different getters. The bridge asks for both and CHECKS
            % THAT THEY AGREE on the state count before pairing Q with pi. An
            % unchecked pairing is the one way this goes wrong in silence: a Q
            % from one enumeration indexed by another's states is not a chain,
            % and every consumer reading a row of the space would read the wrong
            % state rather than get an error.
            gen = CPPLINE.analysisViaCpp(solverName, model, options, 'gen');
            st = CPPLINE.analysisViaCpp(solverName, model, options, 'states');
            n = 0;
            if isfield(gen, 'size')
                n = double(gen.size);
            end
            space = [];
            if isfield(st, 'space')
                space = double(st.space);
            end
            if size(space, 1) ~= n
                line_error(mfilename, sprintf(['line-cli returned a %d-state generator and a ' ...
                    '%d-row state space for the same model; they cannot be paired.'], ...
                    n, size(space, 1)));
            end
            g.infGen = CPPLINE.denseFromTriplets(gen, n);
            g.space = space;
            g.spaceAggr = [];
            if isfield(st, 'spaceAggr')
                g.spaceAggr = double(st.spaceAggr);
            end
            g.pi = [];
            if isfield(st, 'pi')
                g.pi = CPPLINE.jsonNumericList(st.pi)';
            end
            g.eventFilt = {};
            if isfield(gen, 'sync')
                sync = gen.sync;
                if isstruct(sync)
                    sync = num2cell(sync);
                end
                g.eventFilt = cell(1, numel(sync));
                for k = 1:numel(sync)
                    g.eventFilt{k} = CPPLINE.denseFromTriplets(sync{k}, n);
                end
            end
            % Derived START/PREEMPT filtrations, sent as their own keys: they
            % ride on the same arcs as the synchronizations above, so they are
            % kept out of eventFilt (which callers sum as D1). Absent from an
            % older line-cli, in which case the accessors report that this
            % generator carries none rather than reporting zeros as measured.
            g.auxFilt = [];
            if isfield(gen, 'startFilt') || isfield(gen, 'preemptFilt')
                g.auxFilt.start = CPPLINE.auxFiltFromPayload(gen, 'startFilt', n);
                g.auxFilt.preempt = CPPLINE.auxFiltFromPayload(gen, 'preemptFilt', n);
            end
        end

        function F = auxFiltFromPayload(gen, key, n)
            % F = AUXFILTFROMPAYLOAD(GEN, KEY, N)
            % Rebuild an {nstations x nclasses} cell of sparse matrices from the
            % wire's per-(station,class) triplet blocks. Each block carries its
            % own Station and Class fields, 0-based like every other index on
            % this transport.
            F = {};
            if ~isfield(gen, key)
                return
            end
            blocks = gen.(key);
            if isstruct(blocks)
                blocks = num2cell(blocks);
            end
            for b = 1:numel(blocks)
                blk = blocks{b};
                i = double(blk.Station) + 1;
                r = double(blk.Class) + 1;
                F{i,r} = CPPLINE.denseFromTriplets(blk, n); %#ok<AGROW>
            end
            for i = 1:size(F,1)
                for r = 1:size(F,2)
                    if isempty(F{i,r})
                        F{i,r} = sparse(n, n);
                    end
                end
            end
        end

        function M = denseFromTriplets(payload, n)
            % M = DENSEFROMTRIPLETS(PAYLOAD, N)
            % Rebuild an (n x n) sparse matrix from the wire's From/To/Rate
            % triplets, whose indices are 0-based.
            rows = []; cols = []; vals = [];
            if isfield(payload, 'From'),  rows = CPPLINE.jsonNumericList(payload.From); end
            if isfield(payload, 'To'),    cols = CPPLINE.jsonNumericList(payload.To);   end
            if isfield(payload, 'Rate'),  vals = CPPLINE.jsonNumericList(payload.Rate); end
            if numel(rows) ~= numel(cols) || numel(cols) ~= numel(vals)
                line_error(mfilename, sprintf(['line-cli sent %d row indices, %d column ' ...
                    'indices and %d rates for one sparse matrix; the triplets must be of ' ...
                    'one length.'], numel(rows), numel(cols), numel(vals)));
            end
            M = sparse(rows(:) + 1, cols(:) + 1, vals(:), n, n);
        end

        function doc = exportODEs(solverName, model, options, notation)
            % DOC = EXPORTODES(SOLVERNAME, MODEL, OPTIONS, NOTATION)
            % The LaTeX document of the fluid drift, from -a odes. NOTATION is
            % 'scalar' or 'matrix' and REACHES the C++ exporter, so 'matrix'
            % returns the matrix document and not the scalar one under another
            % name.
            if nargin < 4 || isempty(notation)
                notation = 'scalar';
            end
            payload = CPPLINE.analysisViaCpp(solverName, model, options, 'odes', ...
                {'--notation', char(notation)});
            if ~isfield(payload, 'latex')
                line_error(mfilename, 'line-cli''s -a odes payload carries no ''latex'' field.');
            end
            doc = char(payload.latex);
        end

        %% ---- random environment (Environment) -----------------------------

        function [QN, UN, TN, runtime] = getEnvAvg(env, refModel, options)
            % [QN,UN,TN,RUNTIME] = GETENVAVG(ENV, REFMODEL, OPTIONS)
            % Delegate a random-environment solve to line-cli's -s env arm.
            % ENV is the Environment, REFMODEL the stage model whose station and
            % class indexing the returned matrices carry.
            %
            % RN and AN are NOT returned: @SolverENV/getEnsembleAvg defines them
            % as NaN because the environment analyzer computes no response time,
            % and the CLI emits nan for the same reason.
            Tstart = tic;
            binary = CPPLINE.findLineCli();

            tmpDir = tempname;
            if ~mkdir(tmpDir)
                line_error(mfilename, sprintf('could not create the temporary directory ''%s''.', tmpDir));
            end
            cleaner = onCleanup(@() CPPLINE.rmdirQuiet(tmpDir)); %#ok<NASGU>
            modelPath = fullfile(tmpDir, 'model.json');
            linemodel_save(env, modelPath);
            if exist(modelPath, 'file') ~= 2
                line_error(mfilename, 'linemodel_save did not write a model.json for the C++ environment solver.');
            end

            args = {'-f', modelPath, '-i', 'json', '-s', 'env', '-a', 'avg', '-o', 'json'};
            if isstruct(options) && isfield(options, 'method') && ~isempty(options.method) && ...
                    ~strcmp(char(options.method), 'default')
                args = [args, {'--method', char(options.method)}];
            end
            if isstruct(options) && isfield(options, 'arith') && ~isempty(options.arith) && ...
                    ~strcmpi(char(options.arith), 'double')
                args = [args, {'--arith', char(options.arith)}];
            end
            % THE STAGE HORIZON, which @SolverENV passes in options.stagetimespan
            % because it lives on the STAGE SOLVER and not on the ensemble. Left
            % unstated the engine integrates to its own default 100, so the exit
            % averages are a different quadrature -- 6e-3 relative on
            % renv_twostages_repairmen, whose stages ask for [0,1e3].
            % WHICH SOLVER RUNS EACH STAGE. model.json carries the stage
            % Networks and the transition process, and the ensemble's solver
            % choice lives on the SolverENV; the engine defaults to the fluid
            % transient, so an ensemble built on SolverCTMC stages had to be
            % refused outright until this crossed. A cutoff travels with it,
            % since an open stage's chain has to be truncated somewhere.
            if isstruct(options) && isfield(options, 'stagesolver') && ...
                    ~isempty(options.stagesolver)
                args = [args, {'--stage-solver', char(options.stagesolver)}];
                if isfield(options, 'stagecutoff') && ~isempty(options.stagecutoff) && ...
                        isfinite(options.stagecutoff)
                    args = [args, {'--cutoff', sprintf('%d', round(options.stagecutoff))}];
                end
            end
            if isstruct(options) && isfield(options, 'stagetimespan')
                ts = options.stagetimespan;
                if numel(ts) >= 2 && all(isfinite(ts(1:2))) && ts(2) > ts(1)
                    args = [args, {'--tspan', sprintf('%.17g:%.17g', ts(1), ts(2))}];
                end
            end
            results = CPPLINE.runLineCli(binary, args);
            sn = refModel.getStruct();
            [QN, UN, ~, TN] = CPPLINE.avgMatrices(results, sn);
            runtime = toc(Tstart);
        end

        %% ---- LayeredNetwork (LQN) ----------------------------------------

        function [QN, UN, RN, TN, AN, WN, runtime] = getEnsembleAvg(solver, options)
            % [QN,UN,RN,TN,AN,WN,RUNTIME] = GETENSEMBLEAVG(SOLVER, OPTIONS)
            % Delegate a layered getEnsembleAvg() to line-cli, returning the
            % per-LQN-element metric column vectors the MATLAB SolverLN uses.
            %
            % ONE subprocess solves the whole ensemble: letting the MATLAB fixed
            % point run and dispatching each layer separately would spawn one
            % process per layer per iteration.
            Tstart = tic;
            results = CPPLINE.solveLqnViaCpp(solver, options);
            lqn = solver.lqn;
            nidx = numel(lqn.names);
            QN = nan(nidx, 1); UN = nan(nidx, 1); RN = nan(nidx, 1);
            TN = nan(nidx, 1); AN = nan(nidx, 1); WN = nan(nidx, 1);

            if ~isstruct(results) || ~isfield(results, 'rows')
                line_error(mfilename, 'line-cli did not return a layered AvgTable.');
            end
            rows = results.rows;
            if isstruct(rows)
                rows = num2cell(rows);
            end

            % ROWS ARE REORDERED INTO LQN INDEX ORDER, KEYED BY NAME. The two
            % orders are NOT the same: line-cli indexes the elements as the
            % .lqnx declares them, which groups tasks under their processor,
            % while the MATLAB struct numbers tasks in declaration order. On a
            % model whose tasks are declared in one order and hosted in another
            % the two tables carry the same numbers against different rows.
            idxOfName = configureDictionary('string', 'double');
            for idx = 1:nidx
                idxOfName(string(lqn.names{idx})) = idx;
            end
            for k = 1:numel(rows)
                r = rows{k};
                nm = string(r.node);
                if ~isKey(idxOfName, nm)
                    continue
                end
                idx = idxOfName(nm);
                QN(idx) = CPPLINE.jsonScalar(r, 'QLen');
                UN(idx) = CPPLINE.jsonScalar(r, 'Util');
                RN(idx) = CPPLINE.jsonScalar(r, 'RespT');
                WN(idx) = CPPLINE.jsonScalar(r, 'ResidT');
                TN(idx) = CPPLINE.jsonScalar(r, 'Tput');
                % AN stays NaN: the C++ layered result has no arrival-rate
                % vector, and neither does the MATLAB layered table, so the two
                % agree by construction rather than by omission.
            end
            runtime = toc(Tstart);
        end

        function results = solveLqnViaCpp(solver, options)
            % RESULTS = SOLVELQNVIACPP(SOLVER, OPTIONS)
            % Serialize the SolverLN's LayeredNetwork, run line-cli on it and
            % return the parsed object. The wire format is .lqnx unless that
            % would drop a non-reference task's think time, in which case the
            % model.json interchange carries it instead -- see lqnxLossyTasks.
            binary = CPPLINE.findLineCli();
            model = solver.model;
            lossy = CPPLINE.lqnxLossyTasks(solver);
            layerSolver = CPPLINE.lnLayerSolver(solver);

            arith = '';
            if isstruct(options) && isfield(options, 'arith') && ~isempty(options.arith)
                arith = char(options.arith);
            end
            if ~isempty(arith) && ~strcmpi(arith, 'double') && strcmp(layerSolver, 'fluid')
                % line-cli refuses this itself; naming it here says which of the
                % two settings has to change.
                line_error(mfilename, sprintf(['the C++ fluid layer solver is ' ...
                    'double-precision only; arith=''%s'' cannot be combined with Fluid ' ...
                    'layers.'], arith));
            end

            tmpDir = tempname;
            if ~mkdir(tmpDir)
                line_error(mfilename, sprintf('could not create the temporary directory ''%s''.', tmpDir));
            end
            cleaner = onCleanup(@() CPPLINE.rmdirQuiet(tmpDir)); %#ok<NASGU>
            if isempty(lossy)
                modelPath = fullfile(tmpDir, 'model.lqnx');
                model.writeXML(modelPath);
                if exist(modelPath, 'file') ~= 2
                    line_error(mfilename, 'writeXML did not write a model.lqnx for the C++ layered solver.');
                end
                fmtArgs = {'-i', 'lqnx'};
            else
                % No -i: line-cli takes the layered path off the document's own
                % type, and naming -i json would send it to the Network reader.
                modelPath = fullfile(tmpDir, 'model.json');
                linemodel_save(model, modelPath);
                if exist(modelPath, 'file') ~= 2
                    line_error(mfilename, 'linemodel_save did not write a model.json for the C++ layered solver.');
                end
                fmtArgs = {};
            end

            % -s ln and -s ln.mva are the same engine in line-cli (the layer
            % solver comes from --layer-solver, not from the token), so the
            % token stays 'ln' and the layer engine is stated once, explicitly.
            args = [{'-f', modelPath}, fmtArgs, {'-s', 'ln', '-a', 'avg', '-o', 'json', ...
                '--layer-solver', layerSolver}];
            if ~isempty(arith) && ~strcmpi(arith, 'double')
                args = [args, {'--arith', arith}];
            end
            args = [args, CPPLINE.lnKnobs(options)];
            results = CPPLINE.runLineCli(binary, args);
        end

        function layerSolver = lnLayerSolver(solver)
            % LAYERSOLVER = LNLAYERSOLVER(SOLVER)
            % Resolve --layer-solver from the layer solver the MATLAB SolverLN
            % was built with, refusing anything the port has no layer engine
            % for. The layer solver is what the fixed point is a fixed point
            % of, so substituting MVA for a CTMC layer answers a different
            % question. NC and SSA layers ARE carried (solve_layer_nc /
            % solve_layer_ssa in solver_ln.h); refusing them here was a stale
            % claim that left lqn_twotasks and lqn_ofbiz with no table at all.
            layerSolver = 'mva';
            name = '';
            if isprop(solver, 'solverFactory') && ~isempty(solver.solverFactory) && ...
                    isa(solver.solverFactory, 'function_handle')
                probe = solver.solverFactory(CPPLINE.probeNetwork());
                name = class(probe);
            end
            if isempty(name)
                % No factory recorded: the MATLAB default layer engine is MVA.
                return
            end
            short = upper(name);
            if strncmp(short, 'SOLVER', 6)
                short = short(7:end);
            end
            switch short
                case 'MVA'
                    layerSolver = 'mva';
                case {'FLD', 'FLUID'}
                    layerSolver = 'fluid';
                case {'NC', 'COMOM'}
                    layerSolver = 'nc';
                case 'SSA'
                    layerSolver = 'ssa';
                otherwise
                    line_error(mfilename, sprintf(['lang=''cpp'' runs the LQN layers under ' ...
                        'SolverMVA, SolverNC, SolverFluid or SolverSSA; this SolverLN builds ' ...
                        'a ''%s'' layer solver, which the C++ port has no layered engine ' ...
                        'for. Solve it with lang=''matlab'', or build the SolverLN with one ' ...
                        'of those layer factories.'], name));
            end
        end

        function model = probeNetwork()
            % MODEL = PROBENETWORK()
            % A minimal but valid closed network, used only to ask a SolverLN
            % layer factory which solver class it builds. An adaptive factory
            % inspects the layer it is given, so an empty Network would not
            % survive the call.
            model = Network('probe');
            delay = Delay(model, 'ProbeDelay');
            queue = Queue(model, 'ProbeQueue', SchedStrategy.PS);
            jobclass = ClosedClass(model, 'ProbeClass', 1, delay);
            delay.setService(jobclass, Exp(1));
            queue.setService(jobclass, Exp(1));
            model.link(Network.serialRouting(delay, queue));
        end

        function offenders = lqnxLossyTasks(solver)
            % OFFENDERS = LQNXLOSSYTASKS(SOLVER)
            % Name the tasks whose think time an .lqnx serialization of this
            % model would drop, i.e. the non-reference ones that carry one.
            %
            % The lqnx schema accepts think time on a REFERENCE task only, and
            % writeXML reports the loss and writes the file anyway (lqns would
            % reject it otherwise). LINE nonetheless gives a non-reference
            % task's think time to its callers as a delay, so a model carrying
            % one solves to DIFFERENT numbers through that transport: on
            % gallery_lqn_basic, T3's think time caps its throughput at
            % multiplicity/think = 25/4 and drops it from 66.4 to 6.22. A
            % non-empty return therefore does not refuse the solve -- it routes
            % it through linemodel_save, whose reader takes think time on any
            % task.
            model = solver.model;
            offenders = {};
            for t = 1:numel(model.tasks)
                task = model.tasks{t};
                if SchedStrategy.fromText(task.scheduling) == SchedStrategy.REF
                    continue
                end
                if ~isempty(task.thinkTimeMean) && task.thinkTimeMean > 0
                    offenders{end+1} = task.name; %#ok<AGROW>
                end
            end
        end

        function args = lnKnobs(options)
            % ARGS = LNKNOBS(OPTIONS)
            % Translate the SolverLN options into layered-path CLI flags,
            % refusing every setting the flag set cannot carry. A knob with no
            % CLI counterpart is an error and not a silent drop: relax,
            % relax_factor and layering all change the fixed point the MATLAB
            % solver converges to.
            args = {};
            if isempty(options) || ~isstruct(options)
                return
            end
            if isfield(options, 'method') && ~isempty(options.method) && ...
                    ~strcmp(char(options.method), 'default')
                % --method is refused by line-cli on the layered path (it names
                % an algorithm inside a Network solver), and the MATLAB LN
                % methods that are not 'default' are separate engines.
                line_error(mfilename, sprintf(['lang=''cpp'' serves the layered solver''s ' ...
                    'default method only; options.method=''%s'' selects a different layered ' ...
                    'engine (the C++ layered path takes no --method). Use lang=''matlab'' ' ...
                    'for it.'], char(options.method)));
            end
            iterTol = CPPLINE.overridden(options, 'iter_tol', 5e-3);
            if ~isempty(iterTol)
                args = [args, {'--iter_tol', CPPLINE.num2arg(iterTol)}];
            end
            iterMax = CPPLINE.overridden(options, 'iter_max', 200);
            if ~isempty(iterMax) && iterMax > 0
                args = [args, {'--iter_max', sprintf('%d', round(iterMax))}];
            end
            if isfield(options, 'config') && isstruct(options.config)
                cfg = options.config;
                if isfield(cfg, 'interlocking') && ~isempty(cfg.interlocking) && ~cfg.interlocking
                    args = [args, {'--no-interlocking'}];
                end
                if isfield(cfg, 'layering') && ~isempty(cfg.layering) && ...
                        ~strcmpi(char(cfg.layering), 'srvn')
                    line_error(mfilename, sprintf(['lang=''cpp'' implements the ''srvn'' ' ...
                        'layering (one submodel per server); config.layering=''%s'' builds a ' ...
                        'different layer decomposition, which the C++ port does not carry. ' ...
                        'Use lang=''matlab'' for it.'], char(cfg.layering)));
                end
            end
        end

        function v = jsonScalar(s, field)
            % V = JSONSCALAR(S, FIELD)
            % A numeric field of a jsondecode'd row, NaN when absent or null.
            v = NaN;
            if isstruct(s) && isfield(s, field)
                x = s.(field);
                if ischar(x) || isstring(x)
                    v = CPPLINE.jsonNumericScalar(x);
                elseif isnumeric(x) && ~isempty(x)
                    v = double(x(1));
                end
            end
        end

    end
end
