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

% Check global verbosity level to determine output suppression
suppressOutput = isempty(LINEVerbose) || LINEVerbose ~= VerboseLevel.DEBUG;

if suppressOutput
    % Suppress output unless in DEBUG mode
    if ispc
        cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > nul 2>&1'];
    elseif isunix
        cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > /dev/null 2>&1'];
    else
        cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > /dev/null 2>&1'];
    end
else
    % Allow output in DEBUG mode
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'"'];
end

[status] = system(cmd);
if  status > 0
    if suppressOutput
        if ispc
            cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > nul 2>&1'];
        else
            cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > /dev/null 2>&1'];
        end
    else
        cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'"'];
    end
    [status] = system(cmd);
    if status > 0
        rt = java.lang.Runtime.getRuntime();
        rt.exec(cmd);
    end
end
end