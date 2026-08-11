function jar_path = lineDownloadJAR(verbose)
% LINEDOWNLOADJAR Download jline.jar from SourceForge if not found locally.
%
% jar_path = lineDownloadJAR()
% jar_path = lineDownloadJAR(verbose)
%
% Downloads jline.jar from SourceForge (https://line-solver.sourceforge.net/latest/jline.jar)
% and places it in the common/ directory (auto-discovered or current repo).
% Returns the path to the JAR after download.
%
% Args:
%   verbose (optional): If true, print status messages (default: true)
%
% Returns:
%   jar_path: Full path to jline.jar
%
% Usage:
%   jar_path = lineDownloadJAR();  % Auto-download and set up
%   jar_path = lineDownloadJAR(false);  % Silent mode

if nargin < 1 || isempty(verbose)
    verbose = true;
end

% 1. Try to locate existing jline.jar
existing_jar = which('jline.jar');
if ~isempty(existing_jar) && isfile(existing_jar)
    if verbose
        fprintf('jline.jar already found at: %s\n', existing_jar);
    end
    jar_path = existing_jar;
    return;
end

% 2. Try to find/create common/ directory
repo_root = lineGetRepoRoot();
if isempty(repo_root)
    repo_root = pwd;
end

common_dir = fullfile(repo_root, 'common');
if ~isdir(common_dir)
    mkdir(common_dir);
    if verbose
        fprintf('Created directory: %s\n', common_dir);
    end
end

jar_path = fullfile(common_dir, 'jline.jar');

% 3. Check if already downloaded to common/
if isfile(jar_path)
    if verbose
        fprintf('jline.jar found in common/: %s\n', jar_path);
    end
    addPathIfNeeded(common_dir);
    return;
end

% 4. Get version from GlobalConstants.java (works on all OS)
repo_root = lineGetRepoRoot();
gc_file = fullfile(repo_root, 'jar', 'src', 'main', 'java', 'jline', 'GlobalConstants.java');

version = '';
if isfile(gc_file)
    try
        content = fileread(gc_file);
        % Extract version string: Version = "X.Y.Z"
        % Use character arrays for cross-platform compatibility
        tokens = regexp(content, 'Version\s*=\s*"([^"]+)"', 'tokens');
        if ~isempty(tokens) && ~isempty(tokens{1})
            version = char(tokens{1}{1});
        end
    catch
        % Silently continue if version extraction fails
    end
end

% 5. Download versioned jar from SourceForge
if isempty(version)
    % Fallback: download generic jline.jar if version detection fails
    versioned_jar_path = '';
    jar_url = 'https://line-solver.sourceforge.net/latest/jline.jar';
    jar_filename = 'jline.jar';
    if verbose
        fprintf('Version detection failed; will download generic jline.jar\n');
    end
else
    jar_filename = sprintf('jline-%s.jar', version);
    versioned_jar_path = fullfile(common_dir, jar_filename);
    jar_url = sprintf('https://line-solver.sourceforge.net/latest/%s', jar_filename);
end

if verbose
    fprintf('Downloading %s from: %s\n', jar_filename, jar_url);
    fprintf('Target: %s\n', jar_path);
end

try
    % Check if versioned jar already exists (cross-platform)
    if ~isempty(versioned_jar_path) && isfile(versioned_jar_path)
        if verbose
            fprintf('Found existing %s\n', jar_filename);
        end
        if ~isfile(jar_path)
            % Copy versioned jar to generic jline.jar (works on all OS)
            copyfile(versioned_jar_path, jar_path);
            if verbose
                fprintf('Linked: %s -> jline.jar\n', jar_filename);
            end
        end
    else
        % Download from SourceForge (websave works on all OS)
        options = weboptions('Timeout', 60, 'ContentType', 'auto');
        websave(jar_path, jar_url, options);

        if ~isfile(jar_path)
            error('lineDownloadJAR:DownloadFailed', 'Download completed but file not found');
        end

        file_info = dir(jar_path);
        if verbose
            fprintf('Downloaded jline.jar (%.1f MB)\n', file_info.bytes / (1024*1024));
        end

        % If we downloaded versioned jar, also copy to common name
        if ~isempty(versioned_jar_path) && ~isequal(jar_path, versioned_jar_path)
            copyfile(jar_path, versioned_jar_path);
            if verbose
                fprintf('Also saved as: %s\n', jar_filename);
            end
        end
    end

    if ~isfile(jar_path)
        error('lineDownloadJAR:DownloadFailed', 'Failed to set up jline.jar');
    end

    % Add to MATLAB path if not already there
    addPathIfNeeded(common_dir);

catch ME
    % Clean up partial files if download failed (safe on all OS)
    if isfile(jar_path)
        try
            delete(jar_path);
        catch
        end
    end
    if ~isempty(versioned_jar_path) && isfile(versioned_jar_path)
        try
            delete(versioned_jar_path);
        catch
        end
    end

    error('lineDownloadJAR:DownloadFailed', ...
        ['Failed to download jline.jar from SourceForge:\n', ...
         ME.message, '\n\n', ...
         'Options:\n', ...
         '1. Download manually: https://line-solver.sourceforge.net/latest/jline.jar\n', ...
         '2. Build locally: cd jar && mvn clean package -P b\n', ...
         '3. Set JLINE_JAR environment variable to point to existing jar']);
end

end

% -----------------------------------------------------------------------
function repo_root = lineGetRepoRoot()
% Try to find the LINE repository root by searching for key markers.
repo_root = '';

% Start from this file's location
this_file = mfilename('fullpath');
search_dir = fileparts(fileparts(fileparts(this_file)));

% Walk up looking for markers (jar/, matlab/, python/ subdirs)
for depth = 0:5
    if isfolder(fullfile(search_dir, 'jar')) && ...
       isfolder(fullfile(search_dir, 'matlab')) && ...
       isfolder(fullfile(search_dir, 'python'))
        repo_root = search_dir;
        return;
    end
    parent = fileparts(search_dir);
    if strcmp(parent, search_dir)
        break;
    end
    search_dir = parent;
end

end

% -----------------------------------------------------------------------
function addPathIfNeeded(dir_path)
% Add directory to MATLAB path if not already there
current_paths = path;
if ~contains(current_paths, dir_path)
    addpath(dir_path);
end
end
