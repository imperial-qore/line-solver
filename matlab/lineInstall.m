cwd = fileparts(mfilename('fullpath'));
addpath(genpath(cwd));
w=warning('query');
warning on
disp('Checking JAVA...')
[status,~] = system('java -version');
hasWarnings = false;
v = ver;
if status ~= 0
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
[status,~] = system('lqns --help');
if (isunix & status == 127) | (ispc & status > 0) %#ok<AND2> % command not found
    warning('WARNING: LQNS is not installed, so SolverLQNS cannot run. It needs the lqns, lqsim and qnsolver commands. Download them at: https://github.com/layeredqueuing/V6')
    hasWarnings = true;
end
disp('Checking symbolic backend (line-sage-rest)...')
[dstatus,~] = system('docker info');
if dstatus ~= 0
    warning(['WARNING: Docker is not available, so the SageMath symbolic backend cannot start. ', ...
        'It is required by SolverCTMC/SolverFluid symbolic methods (config.symbolic=''sage''). ', ...
        'Install Docker, then run: docker pull imperialqore/line-sage-rest:latest'])
    hasWarnings = true;
else
    [istatus,iresult] = system('docker images -q imperialqore/line-sage-rest');
    if istatus ~= 0 || isempty(strtrim(iresult))
        warning(['WARNING: the line-sage-rest image is not present locally, this may be required by some LINE methods. ', ...
            'Pull it with: docker pull imperialqore/line-sage-rest:latest'])
        hasWarnings = true;
    end
end
disp('Checking JMT...')
lineStart;
% jmtFindPath, not jmtGetPath: asking whether JMT is installed must not be
% what installs it. jmtGetPath downloads 50MB when the jar is absent, so the
% check used to acquire the dependency it was reporting on.
jmtPath = jmtFindPath;
if isempty(jmtPath)
    warning(['WARNING: JMT.jar was not found, so the JMT simulation solver cannot run yet. ', ...
        'It is about 50MB and is downloaded automatically on the first call to the solver.'])
    hasWarnings = true;
else
    disp(['  ', fullfile(jmtPath, 'JMT.jar')])
end
warning(w);
if hasWarnings
    disp('Completed. LINE has warnings.')
else
    disp('Success. LINE is ready to use.')
end
