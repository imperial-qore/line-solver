function cli_path = lineDownloadLineCli(verbose)
% LINEDOWNLOADLINECLI Download the line-cli binary from SourceForge if absent.
%
% cli_path = lineDownloadLineCli()
% cli_path = lineDownloadLineCli(verbose)
%
% Fetches the line-cli build published for this platform under
% https://line-solver.sourceforge.net/latest/ and places it in the common/
% directory of this checkout or release, marked executable. It is the C++ twin
% of lineDownloadJAR, and the twin of find_line_cli's fetch in the native Python
% package (python/line_solver/solvers/cpp_dispatch.py).
%
% WHAT IS PUBLISHED IS A ZIP, and the uncompressed name is its fallback. The
% binary is ~94 MB stripped and ~36 MB deflated, so the archive is what a
% first-time lang='cpp' user waits through; the uncompressed copy stays
% published for LINE releases older than the archive and is what this function
% falls back to if the .zip is missing. MATLAB's unzip does carry the exec bit
% across and Python's zipfile does not; the chmod below is unconditional either
% way, exactly as it already was for the uncompressed copy.
%
% A BINARY IS NOT A JAR. jline.jar runs wherever a JVM does; line-cli is one ELF
% per host triple, so the fetch is gated on a name the platform table below
% knows rather than attempted and diagnosed afterwards from an exec failure.
% upload-line-cli.sh publishes exactly the names this table asks for.
%
% Args:
%   verbose (optional): If true, print status messages (default: true)
%
% Returns:
%   cli_path: Full path to the downloaded line-cli
%
% Errors when this platform has no published build, when LINE_CLI_DOWNLOAD=0
% refuses the download, or when the fetch fails. It never returns '': a caller
% that cannot tell "no binary" from "a binary this call failed to obtain" is
% how a lang='cpp' request comes back with numbers from another engine.

if nargin < 1 || isempty(verbose)
    verbose = true;
end

exeName = 'line-cli';
if ispc
    exeName = 'line-cli.exe';
end

repo_root = lineCliRepoRoot();
common_dir = fullfile(repo_root, 'common');
cli_path = fullfile(common_dir, exeName);

% Already there (a concurrent session may have fetched it since the caller looked).
if isfile(cli_path)
    if verbose
        fprintf('line-cli found in common/: %s\n', cli_path);
    end
    return;
end

tag = lineCliPlatformTag();
if isempty(tag)
    line_error(mfilename, sprintf(['no ''line-cli'' build is published for %s/%s, so ' ...
        'lang=''cpp'' cannot run here. Build one with ''cpp/make.sh -O'' and set ' ...
        'LINE_CLI_BINARY to it, or put it on PATH.'], computer('arch'), computer));
end

if strcmp(strtrim(getenv('LINE_CLI_DOWNLOAD')), '0')
    line_error(mfilename, ['the C++ solver binary ''line-cli'' was not found and ' ...
        'LINE_CLI_DOWNLOAD=0 refuses the download. Set LINE_CLI_BINARY to an existing ' ...
        'binary, put it on PATH, or unset LINE_CLI_DOWNLOAD.']);
end

if ~isfolder(common_dir)
    mkdir(common_dir);
end

% Versioned name first, un-versioned as the fallback, as lineDownloadJAR
% resolves jline.jar: a solver binary that does not match the release driving it
% answers with another engine's feature set, and the un-versioned name serves
% whatever release is current rather than this one.
%
% THE VERSION IS THE PRIMARY KEY AND THE COMPRESSION ONLY THE SECONDARY ONE, so
% the versioned uncompressed copy is tried BEFORE the rolling archive. Ordering
% the two .zip names first would look like a size optimization and would in fact
% silently pin a different release whenever a server carries the rolling archive
% and this release's binary only uncompressed.
version = lineCliVersion(repo_root);
names = {};
if ~isempty(version)
    names{end+1} = sprintf('line-cli-%s-%s.zip', version, tag);
    names{end+1} = sprintf('line-cli-%s-%s', version, tag);
end
names{end+1} = sprintf('line-cli-%s.zip', tag);
names{end+1} = sprintf('line-cli-%s', tag);

% Downloaded under a temporary name and renamed only once complete: a
% half-written binary left in common/ is found by the next call and handed to
% the shell as a solver. An archive is extracted into a temporary directory
% under the same rule, and only the extracted binary is renamed into place.
reasons = {};
for k = 1:numel(names)
    name = names{k};
    isZip = numel(name) > 4 && strcmpi(name(end-3:end), '.zip');
    url = sprintf('https://line-solver.sourceforge.net/latest/%s', name);
    % MATLAB's unzip APPENDS '.zip' to a name that lacks it, so the archive is
    % downloaded under a name that already ends in it rather than under '.part'.
    if isZip
        tmp_path = [cli_path '.part.zip'];
    else
        tmp_path = [cli_path '.part'];
    end
    exdir = '';
    if verbose
        fprintf('Downloading %s from SourceForge...\n', name);
    end
    try
        opts = weboptions('Timeout', 300, 'ContentType', 'binary');
        websave(tmp_path, url, opts);
        info = dir(tmp_path);
        if isempty(info) || info(1).bytes == 0
            error('lineDownloadLineCli:Empty', 'download completed but the file is empty');
        end
        if isZip
            % The directory is created HERE and not inside the helper: MATLAB
            % assigns no output when a function errors, so a helper that made
            % its own would leak a 94 MB extraction tree on every failure.
            exdir = tempname;
            mkdir(exdir);
            src_path = lineCliExtract(tmp_path, exeName, exdir);
        else
            src_path = tmp_path;
        end
        movefile(src_path, cli_path, 'f');
        if ~ispc
            fileattrib(cli_path, '+x');
        end
        if verbose
            final = dir(cli_path);
            if isZip
                fprintf('Downloaded %s (%.1f MB) and extracted %.1f MB to %s\n', ...
                    name, info(1).bytes / 1048576, final(1).bytes / 1048576, cli_path);
            else
                fprintf('Downloaded %s (%.1f MB) to %s\n', name, info(1).bytes / 1048576, cli_path);
            end
        end
        lineCliCleanup(tmp_path, exdir);
        return;
    catch ME
        reasons{end+1} = sprintf('%s: %s', url, ME.message); %#ok<AGROW>
        lineCliCleanup(tmp_path, exdir);
    end
end

line_error(mfilename, sprintf(['the C++ solver binary ''line-cli'' was not found and ' ...
    'the download failed (%s). Set LINE_CLI_BINARY to its path, put it on PATH, or ' ...
    'build it with ''cpp/make.sh -O''.'], strjoin(reasons, '; ')));

end

% -----------------------------------------------------------------------
function bin_path = lineCliExtract(zip_path, exeName, exdir)
% Extract the single `line-cli` entry of a published archive into exdir, which
% the caller owns and removes.
%
% The magic-byte check is not ceremony: a web host that answers a missing file
% with a 200 and an HTML error page hands websave a perfectly good file, and
% unzip's complaint about it reads as a corrupt download rather than as a name
% that is not published. Four bytes distinguish the two.
fid = fopen(zip_path, 'r');
if fid < 0
    error('lineDownloadLineCli:Unreadable', 'the downloaded archive could not be opened');
end
magic = fread(fid, 4, '*uint8')';
fclose(fid);
if ~isequal(magic, uint8([80 75 3 4]))
    error('lineDownloadLineCli:NotAZip', ...
        'the download is not a zip archive (the server may have served an error page)');
end

files = unzip(zip_path, exdir);

bin_path = '';
for k = 1:numel(files)
    [~, base, ext] = fileparts(files{k});
    if strcmpi([base ext], exeName)
        bin_path = files{k};
        break
    end
end
if isempty(bin_path)
    error('lineDownloadLineCli:NoEntry', ...
        'the archive holds no ''%s'' entry (found: %s)', exeName, strjoin(files, ', '));
end
end

% -----------------------------------------------------------------------
function lineCliCleanup(tmp_path, exdir)
% Remove the partial download and any extraction directory, on both paths out of
% the fetch loop. A leftover .part.zip is harmless; a leftover extraction tree is
% a second copy of a 94 MB binary in the temporary directory.
if ~isempty(tmp_path) && isfile(tmp_path)
    try
        delete(tmp_path);
    catch
    end
end
if ~isempty(exdir) && isfolder(exdir)
    try
        rmdir(exdir, 's');
    catch
    end
end
end

% -----------------------------------------------------------------------
function tag = lineCliPlatformTag()
% Host triple of the published line-cli builds, or '' when there is none.
% One entry per build upload-line-cli.sh actually uploads.
tag = '';
if isunix && ~ismac && strcmp(computer('arch'), 'glnxa64')
    tag = 'linux-x86_64';
end
end

% -----------------------------------------------------------------------
function repo_root = lineCliRepoRoot()
% Root of this checkout or release: the parent of the matlab/ tree. It is the
% directory common/ sits in under both layouts, which the jar/ + python/ marker
% walk is not -- a line-<ver>-matlab.zip has neither.
repo_root = fileparts(lineRootFolder());
end

% -----------------------------------------------------------------------
function version = lineCliVersion(repo_root)
% Version of this LINE tree, as a char row vector, or '' when unknown.
% lineStart.m carries LINE_VERSION in every layout; GlobalConstants.java exists
% only in a development checkout.
version = '';
gc_file = fullfile(repo_root, 'jar', 'src', 'main', 'java', 'jline', 'GlobalConstants.java');
start_file = fullfile(repo_root, 'matlab', 'lineStart.m');
candidates = {gc_file, start_file};
patterns = {'Version\s*=\s*"([^"]+)"', 'LINE_VERSION\s*=\s*''([^'']+)'''};
for k = 1:numel(candidates)
    if ~isfile(candidates{k})
        continue
    end
    try
        tokens = regexp(fileread(candidates{k}), patterns{k}, 'tokens');
        if ~isempty(tokens) && ~isempty(tokens{1})
            version = char(tokens{1}{1});
            return
        end
    catch
        % Try the next candidate rather than failing the download outright.
    end
end
end
