function ret = jsimgOpen()
% RET = JSIMGOPEN()
% Open an empty JSIMgraph session

javaCmd = line_java_cmd(mfilename);

cmd = [javaCmd,' -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain'];
[status, cmdout] = system(cmd);
if status > 0
    % Retry on an old JVM that needs the module system relaxed.
    firstOut = cmdout;
    cmdRetry = [javaCmd,' --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain'];
    [status, cmdout] = system(cmdRetry); %#ok<ASGLU>
    if status > 0
        line_error(mfilename, sprintf(['JSIMgraph could not be started (exit code %d).\n', ...
            'Command: %s\nOutput: %s'], status, cmd, strtrim(firstOut)));
    end
end
ret = status;
end
