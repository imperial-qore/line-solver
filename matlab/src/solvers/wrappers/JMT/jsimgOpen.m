function ret = jsimgOpen()

if ispc
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain > nul 2>&1'];
elseif isunix
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain > /dev/null'];
else
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain > /dev/null'];
end
[status] = system(cmd);
if  status > 0
    cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimgraph.mainGui.JSIMGraphMain'];
    [status] = system(cmd);
    if status > 0
        rt = java.lang.Runtime.getRuntime();
        rt.exec(cmd);
    end
end
end

