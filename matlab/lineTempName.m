function tmpname = lineTempName(solvername)
% LINETEMPNAME Generate a temporary file/directory name in the workspace folder
%
% Files are prefixed with 'tmp_' to identify them as temporary files.
% Uses system temp directory to avoid filesystem interference (e.g., Dropbox).

if nargin >= 1
    basedir = fullfile(tempdir, 'line_workspace', solvername);
else
    basedir = fullfile(tempdir, 'line_workspace');
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