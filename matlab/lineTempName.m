function tmpname = lineTempName(solvername, dockerMountable)
% LINETEMPNAME Generate a temporary file/directory name in the workspace folder
%
% Files are prefixed with 'tmp_' to identify them as temporary files.
% Uses system temp directory to avoid filesystem interference (e.g., Dropbox).
%
% LINETEMPNAME(SOLVERNAME, TRUE) roots the workspace under the user's HOME
% instead of the system temp dir. Snap-confined Docker cannot bind-mount the
% system temp dir (typically /tmp), so backends that dispatch through Docker
% (e.g. LQNS) request a HOME-based, mount-accessible location.

baseroot = tempdir;
if nargin >= 2 && ~isempty(dockerMountable) && dockerMountable
    homedir = getenv('HOME');
    if ~isempty(homedir)
        baseroot = fullfile(homedir, '.line');
    end
end

% LINE_WORKSPACE_ROOT relocates every staged model. run-tests.sh sets it when
% wrapping a solver in a container, so the staging dir is one the container can
% bind-mount; nothing in the solvers needs to know a container is involved.
envRoot = getenv('LINE_WORKSPACE_ROOT');
if ~isempty(envRoot)
    baseroot = strtrim(envRoot);
end

if nargin >= 1 && ~isempty(solvername)
    basedir = fullfile(baseroot, 'line_workspace', solvername);
else
    basedir = fullfile(baseroot, 'line_workspace');
end

if ~exist(basedir, 'dir')
    mkdir(basedir);
end

% Generate temp name and add 'tmp_' prefix to the filename part
rawname = tempname(basedir);
[parentdir, filename] = fileparts(rawname);
tmpname = fullfile(parentdir, ['tmp_', filename]);

if ~exist(tmpname, 'dir')
    mkdir(tmpname);
end
end