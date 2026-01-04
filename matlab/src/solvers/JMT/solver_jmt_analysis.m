function [QN,UN,RN,TN,CN,XN,runtime] = solver_jmt_analysis(sn, options)
% [QN,UN,RN,TN,CN,XN,runtime] = SOLVER_JMT_ANALYSIS(SN, OPTIONS)

self.writeJSIM(sn);

cmd = ['java -cp "',getJMTJarPath(self),filesep,'JMT.jar" jmt.commandline.Jmt sim "',getFilePath(self),'jsim',filesep,getFileName(self),'.jsim" -seed ',num2str(options.seed)];
if options.verbose
    line_printf('JMT model: %s\n',[getFilePath(self),'jsim',filesep,getFileName(self),'.jsim']);
    line_printf('JMT command: %s\n',cmd);
end

status = system(cmd);
if  status > 0
    cmd = ['java -cp "',getJMTJarPath(self),filesep,'JMT.jar" jmt.commandline.Jmt sim "',getFilePath(self),'jsim',filesep,getFileName(self),'.jsim" -seed ',num2str(options.seed),' --illegal-access=permit'];
    [status] = system(cmd);
    if status > 0
        rt = java.lang.Runtime.getRuntime();
        rt.exec(cmd);
    end
end

runtime = toc(Tstart);
self.getResults;
end