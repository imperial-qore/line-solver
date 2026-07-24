function data = solveCli(self, options, extraFlags)
% DATA = SOLVECLI(OPTIONS, EXTRAFLAGS)
%
% Fully JSON-mediated LDES solve. Serializes the model natively to model.json
% (linemodel_save, no in-process Java object marshalling), runs the LDES engine
% as a subprocess exchanging only JSON ("solve model.json -o result.json"), and
% returns the jsondecode'd result struct. Returns [] on any failure.
%
% The engine is either the native GraalVM binary (common/ldes) or, failing that,
% "java -jar common/ldes.jar" (same shaded engine). No JPype/JLINE path remains.
%
% EXTRAFLAGS is an optional cellstr of additional CLI tokens, e.g.
% {'--timespan','0,100','--trajectory'} or {'--export-histogram'}.

if nargin < 3 || isempty(extraFlags)
    extraFlags = {};
end
data = [];

tmpDir = tempname;
if ~mkdir(tmpDir)
    return;
end
cleaner = onCleanup(@() rmdirQuiet(tmpDir)); %#ok<NASGU>
modelPath = fullfile(tmpDir, 'model.json');
resultPath = fullfile(tmpDir, 'result.json');

% Serialize the model to the CLI model.json format natively (no LINE2JLINE).
try
    linemodel_save(self.model, modelPath);
catch ME
    line_debug('LDES solveCli: linemodel_save failed: %s', ME.message);
    return;
end
if exist(modelPath, 'file') ~= 2
    return;
end

% Remote engine, if options.rest_url names an LDES REST server. The container
% owns the engine in that case, so no local binary or jar is required.
restUrl = '';
if isfield(options, 'rest_url') && ~isempty(options.rest_url)
    restUrl = char(options.rest_url);
end

% Resolve the ordered engine command prefixes (native binary, then jar).
runners = {};
if isempty(restUrl)
    runners = SolverLDES.getLdesRunners();
    if isempty(runners)
        line_error(mfilename, ['No runnable LDES engine found (neither the native ' ...
            'common/ldes binary nor a Java runtime for common/ldes.jar).']);
    end
end

% see _kb/06-solver-catalog.md (Wrappers: LDES runner ordering vs stale AOT native binary)
if numel(runners) > 1
    preferJar = false;
    % A flag the AOT native binary predates must run on the jar, which is
    % rebuilt with the sources. --respt-samples is the current instance.
    for fi = 1:numel(extraFlags)
        if ischar(extraFlags{fi}) && strcmp(extraFlags{fi}, '--respt-samples')
            preferJar = true;
            break;
        end
    end
    for ni = 1:numel(self.model.nodes)
        nd = self.model.nodes{ni};
        if isa(nd, 'Source') && ~isempty(nd.markedClasses)
            preferJar = true;
            break;
        end
        if isa(nd, 'Station') && ~isempty(nd.lcdScaling)
            preferJar = true;
            break;
        end
    end
    if preferJar
        runners = runners(end:-1:1);
    end
end

% Flags shared across runners (mirror the Python-native flag set).
flags = sprintf('-s %d --seed %d', round(options.samples), round(options.seed));
if isfield(options, 'method') && (ischar(options.method) || isstring(options.method)) ...
        && ~strcmpi(char(options.method), 'default')
    flags = sprintf('%s --method %s', flags, char(options.method));
end
% Warm-start placement (station-major vector) via --initsol.
if isfield(options, 'init_sol') && ~isempty(options.init_sol)
    v = options.init_sol(:).';
    flags = sprintf('%s --initsol %s', flags, strjoin(arrayfun(@(x) num2str(x, '%.10g'), v, ...
        'UniformOutput', false), ','));
end
% Discrete-time (slotted) mode. --slotlength implies --slotted on the CLI side,
% but both are emitted when the length is non-default so the command line states
% the intent explicitly.
slottedOn = logical(ldesOpt(options, 'slotted', false));
if slottedOn
    flags = sprintf('%s --slotted', flags);
    slotLen = ldesOpt(options, 'slotlength', 1);
    if slotLen ~= 1
        flags = sprintf('%s --slotlength %.10g', flags, slotLen);
    end
end
% Transient (warmup) filter and confidence-interval estimator. Only non-default
% values are emitted, so a default run keeps a minimal command line that an
% older AOT native binary still understands.
tranfilter = ldesOpt(options, 'tranfilter', 'mser5');
if ischar(tranfilter) || isstring(tranfilter)
    tranfilter = char(tranfilter);
    if ~strcmpi(tranfilter, 'mser5')
        flags = sprintf('%s --tranfilter %s', flags, tranfilter);
    end
end
warmupfrac = ldesOpt(options, 'warmupfrac', 0.2);
if isnumeric(warmupfrac) && isscalar(warmupfrac) && isfinite(warmupfrac) && warmupfrac ~= 0.2
    flags = sprintf('%s --warmupfrac %.10g', flags, warmupfrac);
end
mserbatch = ldesOpt(options, 'mserbatch', 5);
if isnumeric(mserbatch) && isscalar(mserbatch) && mserbatch ~= 5
    flags = sprintf('%s --mserbatch %d', flags, round(mserbatch));
end
cimethod = ldesOpt(options, 'cimethod', 'obm');
if ischar(cimethod) || isstring(cimethod)
    cimethod = char(cimethod);
    if ~strcmpi(cimethod, 'obm')
        flags = sprintf('%s --cimethod %s', flags, cimethod);
    end
end
obmoverlap = ldesOpt(options, 'obmoverlap', 0.5);
if isnumeric(obmoverlap) && isscalar(obmoverlap) && isfinite(obmoverlap) && obmoverlap ~= 0.5
    flags = sprintf('%s --obmoverlap %.10g', flags, obmoverlap);
end
ciminbatch = ldesOpt(options, 'ciminbatch', 10);
if isnumeric(ciminbatch) && isscalar(ciminbatch) && ciminbatch ~= 10
    flags = sprintf('%s --ciminbatch %d', flags, round(ciminbatch));
end
ciminobs = ldesOpt(options, 'ciminobs', 100);
if isnumeric(ciminobs) && isscalar(ciminobs) && ciminobs ~= 100
    flags = sprintf('%s --ciminobs %d', flags, round(ciminobs));
end
spectralfrac = ldesOpt(options, 'spectrallowfreqfrac', 0.25);
if isnumeric(spectralfrac) && isscalar(spectralfrac) && isfinite(spectralfrac) && spectralfrac ~= 0.25
    flags = sprintf('%s --spectrallowfreqfrac %.10g', flags, spectralfrac);
end
% Convergence-based stopping. The tolerance is options.config.cnvgtol, the field
% SolverOptions declares for LDES; it is not options.iter_tol, which is the
% fixed-point tolerance of the analytical solvers and has an unrelated default.
if logical(ldesOpt(options, 'cnvgon', false))
    flags = sprintf('%s --cnvgon', flags);
    cnvgtol = ldesOpt(options, 'cnvgtol', 0.05);
    if isnumeric(cnvgtol) && isscalar(cnvgtol) && isfinite(cnvgtol) && cnvgtol ~= 0.05
        flags = sprintf('%s --cnvgtol %.10g', flags, cnvgtol);
    end
    cnvgbatch = ldesOpt(options, 'cnvgbatch', 20);
    if isnumeric(cnvgbatch) && isscalar(cnvgbatch) && cnvgbatch ~= 20
        flags = sprintf('%s --cnvgbatch %d', flags, round(cnvgbatch));
    end
    cnvgchk = ldesOpt(options, 'cnvgchk', 0);
    if isnumeric(cnvgchk) && isscalar(cnvgchk) && cnvgchk ~= 0
        flags = sprintf('%s --cnvgchk %d', flags, round(cnvgchk));
    end
end
% Independent replications and the worker pool that runs them. Emitted here for
% every analysis, steady-state included; runTransientJson no longer adds them.
% A caller that already supplied the flag through extraFlags wins.
reps = ldesOpt(options, 'replications', 1);
if isnumeric(reps) && isscalar(reps) && reps > 1 && ~hasFlag(extraFlags, '--replications')
    flags = sprintf('%s --replications %d', flags, round(reps));
    nthreads = ldesOpt(options, 'numthreads', []);
    if isnumeric(nthreads) && isscalar(nthreads) && nthreads > 0 ...
            && ~hasFlag(extraFlags, '--numthreads')
        flags = sprintf('%s --numthreads %d', flags, round(nthreads));
    end
end
% Cooperative wall-clock budget. The SSJ event loop checks it and stops early
% with stoppingReason='max_time'; the subprocess timeout remains the hard bound.
timeout = ldesOpt(options, 'timeout', Inf);
if isnumeric(timeout) && isscalar(timeout) && isfinite(timeout) && timeout > 0 ...
        && ~hasFlag(extraFlags, '--maxtime')
    flags = sprintf('%s --maxtime %.10g', flags, timeout);
end
for i = 1:numel(extraFlags)
    flags = sprintf('%s %s', flags, extraFlags{i});
end

% Remote engine: POST the same model.json and flags the CLI would have used.
if ~isempty(restUrl)
    data = solveRest(restUrl, modelPath, flags, options);
    return;
end

% Try each runner in order; the native binary may lack some reflective features
% (e.g. fork-join MMT serialization), so fall through to the JVM jar on failure.
lastOut = '';
for ri = 1:numel(runners)
    if exist(resultPath, 'file') == 2
        delete(resultPath);
    end
    cmd = sprintf('%s solve "%s" -o "%s" %s', runners{ri}, modelPath, resultPath, flags);
    [status, cmdout] = system(cmd);
    lastOut = cmdout;
    if status == 0 && exist(resultPath, 'file') == 2
        try
            data = jsondecode(fileread(resultPath));
        catch ME
            line_error(mfilename, sprintf('LDES result JSON parse failed: %s', ME.message));
        end
        if isstruct(data) && isfield(data, 'error')
            line_error(mfilename, sprintf('LDES engine error: %s', data.error));
        end
        return;
    end
    line_debug('LDES runner %d/%d failed (status=%d); trying next.', ri, numel(runners), status);
end

line_error(mfilename, sprintf('LDES engine failed on all %d runner(s): %s', ...
    numel(runners), lastOut));
end

function data = solveRest(restUrl, modelPath, flags, options)
% DATA = SOLVEREST(RESTURL, MODELPATH, FLAGS, OPTIONS)
%
% Solve through an LDES REST server (the imperialqore/ldes container). The wire
% format is the same model.json the CLI reads and the same ldes-result document
% it writes, so the returned struct is identical to the subprocess path for a
% fixed seed and sample count.

url = regexprep(char(restUrl), '/+$', '');
if isempty(regexp(url, '/api/v\d+/solve$', 'once'))
    url = [url '/api/v1/solve'];
end

req = struct();
req.model = struct('content', fileread(modelPath), 'base64', false);
% The server accepts the CLI long-form flags verbatim; values never contain
% whitespace, so a whitespace split reproduces the argument vector exactly.
tokens = strsplit(strtrim(flags));
tokens = tokens(~cellfun(@isempty, tokens));
req.flags = tokens;

timeout = 3600;
if isfield(options, 'timeout') && ~isempty(options.timeout) && isfinite(options.timeout)
    timeout = options.timeout;
end
wopts = weboptions('MediaType', 'application/json', 'ContentType', 'json', ...
    'RequestMethod', 'post', 'Timeout', timeout);

try
    resp = webwrite(url, req, wopts);
catch ME
    line_error(mfilename, sprintf('LDES REST request to %s failed: %s', url, ME.message));
end

if ~isstruct(resp) || ~isfield(resp, 'status')
    line_error(mfilename, sprintf('LDES REST server at %s returned an unexpected payload.', url));
end
if ~strcmpi(resp.status, 'ok')
    msg = 'unspecified error';
    if isfield(resp, 'message') && ~isempty(resp.message)
        msg = resp.message;
    end
    detail = '';
    if isfield(resp, 'stderr') && ~isempty(resp.stderr)
        detail = sprintf(' Engine stderr: %s', strtrim(resp.stderr));
    end
    line_error(mfilename, sprintf('LDES REST solve failed: %s.%s', msg, detail));
end

data = resp.result;
if isstruct(data) && isfield(data, 'error')
    line_error(mfilename, sprintf('LDES engine error: %s', data.error));
end
end

function rmdirQuiet(d)
% RMDIRQUIET Remove directory d and contents, ignoring errors.
try
    if exist(d, 'dir')
        rmdir(d, 's');
    end
catch
end
end

function v = ldesOpt(options, name, default)
% V = LDESOPT(OPTIONS, NAME, DEFAULT)
% Read an LDES option under either spelling. SolverOptions declares the engine
% knobs under options.config, but a caller may set the flat field, and
% Solver.parseOptions replaces the defaults wholesale when the caller passes its
% own struct, so neither spelling can be assumed present.
v = default;
if isfield(options, name) && ~isempty(options.(name))
    v = options.(name);
elseif isfield(options, 'config') && isstruct(options.config) ...
        && isfield(options.config, name) && ~isempty(options.config.(name))
    v = options.config.(name);
end
end

function tf = hasFlag(extraFlags, flag)
% TF = HASFLAG(EXTRAFLAGS, FLAG)
% True if the caller already supplied FLAG among the extra CLI tokens.
tf = false;
for i = 1:numel(extraFlags)
    if (ischar(extraFlags{i}) || isstring(extraFlags{i})) ...
            && strcmp(char(extraFlags{i}), flag)
        tf = true;
        return;
    end
end
end

