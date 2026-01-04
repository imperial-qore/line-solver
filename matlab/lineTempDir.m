function tmpdir = lineTempDir
% LINETEMPDIR Generate a temporary directory in the workspace folder
%
% Directories are prefixed with 'tmp_' to identify them as temporary.

basedir = [lineRootFolder, filesep, 'workspace'];

% Generate temp name and add 'tmp_' prefix to the directory name
rawname = tempname(basedir);
[parentdir, dirname] = fileparts(rawname);
tmpdir = [fullfile(parentdir, ['tmp_', dirname]), filesep];

if ~exist(tmpdir, 'dir')
    mkdir(tmpdir);
end
end