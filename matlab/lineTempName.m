function tmpname = lineTempName(solvername)
% LINETEMPNAME Generate a temporary file/directory name in the workspace folder
%
% Files are prefixed with 'tmp_' to identify them as temporary files.

if nargin >= 1
    basedir = [lineRootFolder, filesep, 'workspace', filesep, solvername];
else
    basedir = [lineRootFolder, filesep, 'workspace'];
end

% Generate temp name and add 'tmp_' prefix to the filename part
rawname = tempname(basedir);
[parentdir, filename] = fileparts(rawname);
tmpname = fullfile(parentdir, ['tmp_', filename]);

if ~exist(tmpname, 'dir')
    mkdir(tmpname);
end
end