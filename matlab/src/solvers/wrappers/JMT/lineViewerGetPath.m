function viewerJar = lineViewerGetPath
% VIEWERJAR = LINEVIEWERGETPATH

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get path to line-viewer.jar in common folder
% This function is at dev/src/solvers/wrappers/JMT/lineViewerGetPath.m
% Navigate up to the line-dev.git root directory
current_dir = fileparts(mfilename('fullpath'));           % dev/src/solvers/wrappers/JMT
root_dir = fileparts(fileparts(fileparts(fileparts(fileparts(current_dir))))); % line-dev.git
common_dir = fullfile(root_dir, 'common');
viewerJar = fullfile(common_dir, 'line-viewer.jar');

if ~exist(viewerJar, 'file')
    line_error(mfilename, sprintf('line-viewer.jar not found at %s. Please build it and place it in the common folder.', viewerJar));
end
end
