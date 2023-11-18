function jsimwView(filename)
% JSIMWVIEW(FILENAME)
% Open model in JSIMwiz

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
[path] = fileparts(filename);
if isempty(path)
    filename=[pwd,filesep,filename];
end
if ispc
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > nul 2>&1'];
elseif isunix
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > /dev/null'];
else
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" > /dev/null'];
end
[status] = system(cmd);
if  status > 0
    cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.commandline.Jmt jsimw "',filename,'" '];
    [status] = system(cmd);
    if status > 0
        rt = java.lang.Runtime.getRuntime();
        rt.exec(cmd);
    end
end
end