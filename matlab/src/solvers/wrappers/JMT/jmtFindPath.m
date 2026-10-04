function jmtPath = jmtFindPath
% JMTPATH = JMTFINDPATH  Locate an existing JMT.jar, downloading nothing.
%
% The search half of JMTGETPATH, split out so that a caller that only wants
% to KNOW whether JMT is installed -- the environment check LINEINSTALL --
% does not trigger the 50MB fetch as a side effect of asking. A check must
% not change the environment it reports on.
%
% Returns the folder holding JMT.jar, or '' when there is none.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

current_dir = fileparts(mfilename('fullpath'));           % dev/src/solvers/wrappers/JMT
root_dir = fileparts(fileparts(fileparts(fileparts(fileparts(current_dir))))); % line-dev.git
common_dir = fullfile(root_dir, 'common');
if exist(fullfile(common_dir, 'JMT.jar'), 'file')
    jmtPath = common_dir;
else
    jmtPath = '';
end
end
