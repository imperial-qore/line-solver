function ret = jsimwOpen()

if ispc
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimwiz.JSIMWizMain > nul 2>&1'];
elseif isunix
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimwiz.JSIMWizMain > /dev/null'];
else
    cmd = ['java -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimwiz.JSIMWizMain > /dev/null'];
end
[status] = system(cmd);
if  status > 0
    cmd = ['java --illegal-access=permit -cp "',jmtGetPath,filesep,'JMT.jar" jmt.gui.jsimwiz.JSIMWizMain'];
    [status] = system(cmd);
    if status > 0
        rt = java.lang.Runtime.getRuntime();
        rt.exec(cmd);
    end
end
end

