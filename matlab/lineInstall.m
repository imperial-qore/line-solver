cwd = fileparts(mfilename('fullpath'));
addpath(genpath(cwd));
w=warning('query');
warning on
disp('Checking JAVA...')
[status,~] = system('java');
hasWarnings = false;
v = ver;
if status == 0
    error('ERROR: the Java Runtime Environment (JRE) is not installed, this is required for LINE.')
end
disp('Checking MATLAB toolboxes...')
if ~any(strcmp('Statistics and Machine Learning Toolbox', {v.Name}))
    warning('ERROR: the Statistics and Machine Learning toolbox is not installed, this is required for LINE.')
end
if ~any(strcmp('Optimization Toolbox', {v.Name}))
    warning('ERROR: the Optimization Toolbox is not installed, this is required for LINE.')
end
if ~any(strcmp('Global Optimization Toolbox', {v.Name}))
    warning('ERROR: The Global Optimization Toolbox is not installed, this is required for LINE.')
    hasWarnings = true;
end
if ~any(strcmp('Parallel Computing Toolbox', {v.Name}))
    warning('ERROR: The Parallel Computing Toolbox is not installed, this is required for LINE.')
    hasWarnings = true;
end
if ~any(strcmp('Symbolic Math Toolbox', {v.Name}))
    warning('ERROR: The Symbolic Math Toolbox is not installed, this may be required by some LINE methods.')
    hasWarnings = true;
end
disp('Checking LQNS...')
[status,result] = system('lqns --h');
if (isunix & result == 127) | (ispc & result > 0) %#ok<AND2> % command not found
    warning('WARNING: LQNS is not installed, this may be required by some LINE methods. Download it at: https://github.com/layeredqueuing/dist')
    hasWarnings = true;
end
disp('Checking JMT...')
lineStart;
jmtGetPath;
warning(w);
if hasWarnings
    disp('Completed. LINE has warnings.')
else
    disp('Success. LINE is ready to use.')
end
