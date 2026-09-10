function [status, cmdout] = jmtRun(self, mode, fname, seed, options)
% [STATUS, CMDOUT] = JMTRUN(MODE, FNAME, SEED, OPTIONS)
%
% Run one batch analysis of jmt.commandline.Jmt on the model FNAME and leave
% the result file where the JMT CLI itself would leave it, that is at
% [FNAME,'-result.jsim'] for MODE='sim' and [FNAME,'-result.jmva'] for
% MODE='mva'. Every backend below satisfies that contract, so getResultsJSIM
% and getResultsJMVA parse the same file whichever one ran.
%
% Backend selection, in order:
%   1. options.rest_url non-empty: POST to a JMT REST server (the
%      imperialqore/jmt-rest container). Nothing is executed locally.
%   2. a local JVM plus common/JMT.jar: the default, unchanged. The launcher
%      is resolved by line_java_exe, so a host whose PATH has no java (typical
%      on Windows) still runs through the JRE bundled with MATLAB.
%   3. no local JVM, but Docker is usable: ask the user once per session
%      whether to pull and use the JMT image, and dispatch through it.
%
% see _kb/06-solver-catalog.md (Wrappers: three ways to reach an external binary)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

restUrl = '';
if isfield(options, 'rest_url') && ~isempty(options.rest_url)
    restUrl = char(options.rest_url);
end

if ~isempty(restUrl)
    [status, cmdout] = jmtSolveRest(self, restUrl, mode, fname, seed, options);
    return;
end

javaExe = line_java_exe();
if ~isempty(javaExe)
    cmd = ['"', javaExe, '" -cp "', getJMTJarPath(self), filesep, 'JMT.jar" jmt.commandline.Jmt ', ...
        mode, ' "', fname, '" -seed ', num2str(seed), ' --illegal-access=permit'];
    if options.verbose == VerboseLevel.DEBUG
        line_printf('JMT command: %s\n', cmd);
    end
    [status, cmdout] = system(cmd);
    return;
end

image = jmtDockerImage(options);
if isempty(image)
    line_error(mfilename, ['SolverJMT requires a Java runtime and JMT.jar. No JVM was found ', ...
        'on the path, and Docker is not usable either. Install Java, or start a JMT REST ', ...
        'server and set options.rest_url.']);
end

[status, cmdout] = jmtRunDocker(image, mode, fname, seed, options);
end


function image = jmtDockerImage(options)
% IMAGE = JMTDOCKERIMAGE(OPTIONS)
% Resolve the JMT Docker image to dispatch through, asking the user for
% consent before anything is pulled. Returns '' when Docker is unusable, when
% no image can be obtained, or when the user declines.
%
% The consent is remembered for the session: a sweep over many models must not
% ask once per model.
persistent decision
image = '';

% Docker bind-mount dispatch is supported on unix hosts only, as for LQNS.
if ispc
    return;
end
if unix('docker info >/dev/null 2>&1') ~= 0
    return;
end

if isfield(options, 'config') && isstruct(options.config) ...
        && isfield(options.config, 'container') && ~isempty(options.config.container)
    candidates = {char(options.config.container)};
else
    % LINE_JMT_IMAGE overrides the default image (parity with the Java/Python rows).
    envImage = getenv('LINE_JMT_IMAGE');
    if ~isempty(envImage)
        candidates = {strtrim(envImage)};
    else
        candidates = {'imperialqore/jmt-rest:latest', 'imperialqore/jmt-rest'};
    end
end

% An image already on the host needs no pull and no question.
for i = 1:numel(candidates)
    [st, out] = unix(['docker images -q ', candidates{i}, ' 2>/dev/null']);
    if st == 0 && ~isempty(strtrim(out))
        image = candidates{i};
        return;
    end
end

if strcmp(decision, 'no')
    return;
elseif strcmp(decision, 'yes')
    % Consent given earlier in the session but the pull left nothing behind.
    decision = '';
end

target = candidates{1};
if ~jmtAskDockerConsent(target)
    decision = 'no';
    return;
end
decision = 'yes';

% Storage check before pulling (shared guard used by the LQNS/QNS wrappers).
if ~lineDockerHasStorageFor(target)
    line_warning(mfilename, ['Skipping docker pull of %s: insufficient free space at the ' ...
        'Docker storage location. Install Java, or set options.rest_url to a running ' ...
        'JMT REST server.\n'], target);
    decision = 'no';
    return;
end

line_printf('Pulling %s (this happens once)...\n', target);
if unix(['docker pull ', target]) ~= 0
    line_warning(mfilename, 'Could not pull %s. Install Java, or set options.rest_url to a running JMT REST server.\n', target);
    decision = 'no';
    return;
end
image = target;
end


function tf = jmtAskDockerConsent(image)
% TF = JMTASKDOCKERCONSENT(IMAGE)
% Ask whether the JMT Docker image may be pulled and used. LINE_JMT_DOCKER
% answers for the user in unattended runs: '1' consents, '0' refuses.
env = getenv('LINE_JMT_DOCKER');
if ~isempty(env)
    tf = any(strcmpi(strtrim(env), {'1', 'true', 'yes', 'y'}));
    return;
end

% Headless/unattended runs must never block on a prompt: refuse without asking
% (mirrors the System.console()==null / stdin.isatty() guards in the Java and
% Python rows). LINE_JMT_DOCKER above is how such runs opt in.
if jmtIsHeadless()
    tf = false;
    return;
end

line_printf(['\nSolverJMT needs Java, which was not found on this host.\n', ...
    'Docker is available and can run JMT from the image %s instead.\n'], image);
answer = input('Pull and use that image? [y/N]: ', 's');
% An empty answer, which is also what a non-interactive session returns, is a
% refusal: a solver must never block or download without being asked.
tf = ~isempty(answer) && any(strcmpi(strtrim(answer), {'y', 'yes'}));
end


function tf = jmtIsHeadless()
% TF = JMTISHEADLESS()
% True in an unattended MATLAB session (started with -batch), where input()
% must not be reached. Guarded so it also works on releases without
% batchStartupOptionUsed.
tf = false;
try
    tf = batchStartupOptionUsed;
catch
    % Older MATLAB: no batchStartupOptionUsed; fall back to false and rely on
    % input() returning '' for the non-interactive case.
end
end


function [status, cmdout] = jmtRunDocker(image, mode, fname, seed, options)
% [STATUS, CMDOUT] = JMTRUNDOCKER(IMAGE, MODE, FNAME, SEED, OPTIONS)
% Run the analysis in the JMT container and copy the result back next to
% FNAME, so the caller sees the same layout the local JVM would have produced.
%
% The model is staged under HOME rather than run in place: snap-confined
% Docker cannot bind-mount the system temp dir where the JMT model normally
% lives. In mva mode JMT rewrites the model file itself, which is why the
% staged copy, not the original, is what the container touches.
workdir = lineTempName('jmt-docker', true);
if ~exist(workdir, 'dir')
    mkdir(workdir);
end
cleanup = onCleanup(@() rmdirQuiet(workdir));

[~, base, ext] = fileparts(fname);
staged = fullfile(workdir, [base, ext]);
copyfile(fname, staged);
fileattrib(staged, '+w');

[~, hostUid] = unix('id -u');
[~, hostGid] = unix('id -g');
hostUid = strtrim(hostUid);
hostGid = strtrim(hostGid);
cmd = ['docker run --rm --user ', hostUid, ':', hostGid, ...
    ' -v ', workdir, ':', workdir, ' -w ', workdir, ' ', image, ' ', ...
    mode, ' "', [base, ext], '" -seed ', num2str(seed)];

if options.verbose == VerboseLevel.DEBUG
    line_printf('JMT command: %s\n', cmd);
end
[status, cmdout] = system(cmd);

resultName = [base, ext, '-result.', jmtResultExt(mode)];
staged_result = fullfile(workdir, resultName);
if exist(staged_result, 'file')
    copyfile(staged_result, [fname, '-result.', jmtResultExt(mode)]);
end
end


function ext = jmtResultExt(mode)
% EXT = JMTRESULTEXT(MODE)
% The suffix jmt.commandline.Jmt appends to the model path per analysis.
switch mode
    case 'sim'
        ext = 'jsim';
    case 'mva'
        ext = 'jmva';
    otherwise
        line_error(mfilename, sprintf('Unknown JMT analysis mode: %s', mode));
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
