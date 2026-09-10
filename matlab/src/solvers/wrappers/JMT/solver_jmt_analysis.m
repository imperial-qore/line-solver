function [QN,UN,RN,TN,CN,XN,runtime] = solver_jmt_analysis(sn, options)
% [QN,UN,RN,TN,CN,XN,runtime] = SOLVER_JMT_ANALYSIS(SN, OPTIONS)

self.writeJSIM(sn);

javaCmd = line_java_cmd(mfilename);
cmd = [javaCmd,' -cp "',getJMTJarPath(self),filesep,'JMT.jar" jmt.commandline.Jmt sim "',getFilePath(self),'jsim',filesep,getFileName(self),'.jsim" -seed ',num2str(options.seed)];
if options.verbose
    line_printf('JMT model: %s\n',[getFilePath(self),'jsim',filesep,getFileName(self),'.jsim']);
    line_printf('JMT command: %s\n',cmd);
end

status = system(cmd);
if  status > 0
    cmd = [javaCmd,' -cp "',getJMTJarPath(self),filesep,'JMT.jar" jmt.commandline.Jmt sim "',getFilePath(self),'jsim',filesep,getFileName(self),'.jsim" -seed ',num2str(options.seed),' --illegal-access=permit'];
    [status, cmdout] = system(cmd);
    if status > 0
        line_error(mfilename, sprintf(['JMT simulation failed (exit code %d).\n', ...
            'Command: %s\nOutput: %s'], status, cmd, strtrim(cmdout)));
    end
end

runtime = toc(Tstart);
self.getResults;
end