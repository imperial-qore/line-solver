function jsimwView(filename)
% JSIMWVIEW(FILENAME)
% Open model in JSIMwiz

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

global LINEVerbose;

[path] = fileparts(filename);
if isempty(path)
    filename=[pwd,filesep,filename];
end

% Resolved once here: a missing JVM must be reported as such, not as a shell
% "not recognized" status that the retry below would mistake for a JMT error.
javaCmd = line_java_cmd(mfilename);
verbose = ~isempty(LINEVerbose) && LINEVerbose == VerboseLevel.DEBUG;

% Both attempts capture their output rather than redirecting it away: system
% with two outputs echoes nothing, and the output is what makes a failure
% diagnosable.
cmd = [javaCmd,' -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'"'];
if verbose
    line_printf('JMT command: %s\n',cmd);
end
[status, cmdout] = system(cmd);
if verbose
    line_printf('%s\n',cmdout);
end

if status > 0
    % Retry on an old JVM that needs the module system relaxed. A JVM that does
    % not know the flag rejects it, so the first failure is the one reported.
    firstOut = cmdout;
    cmdRetry = [javaCmd,' --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'"'];
    [status, cmdout] = system(cmdRetry); %#ok<ASGLU>
    if status > 0
        line_error(mfilename, sprintf(['JSIMwiz could not be started (exit code %d).\n', ...
            'Command: %s\nOutput: %s'], status, cmd, strtrim(firstOut)));
    end
end
end
