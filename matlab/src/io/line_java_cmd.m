function javaCmd = line_java_cmd(caller)
% JAVACMD = LINE_JAVA_CMD(CALLER)
%
% Quoted Java launcher to prefix a JVM command line with, or a clear error
% naming the missing runtime when none is reachable. CALLER is the function
% reported in that error, usually mfilename.
%
% Callers must go through here rather than hardcoding 'java': on a host with
% no Java on PATH, the bare string produces "'java' is not recognized" from
% the shell and, if the failure is then retried through
% java.lang.Runtime.exec, an unhandled java.io.IOException with
% CreateProcess error=2 that says nothing about what is actually missing.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1
    caller = 'line_java_cmd';
end

javaExe = line_java_exe();
if isempty(javaExe)
    line_error(caller, ['No Java runtime was found, so JMT cannot be started. ', ...
        'Install a JRE or JDK (https://adoptium.net) and make sure "java" is on the ', ...
        'system PATH, or point LINE at an existing one by setting the JAVA_HOME ', ...
        'environment variable to its installation directory, or LINE_JAVA to the full ', ...
        'path of the java executable (java.exe on Windows). Searched: LINE_JAVA, ', ...
        'JAVA_HOME, MATLAB_JAVA, the JRE bundled with MATLAB, and the PATH.']);
end

% Quoting is mandatory: the bundled JRE lives under a path with spaces on a
% default Windows install.
javaCmd = ['"', javaExe, '"'];
end
