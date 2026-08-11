function modelView(self)
% MODELVIEW() - Open model in JSIMgraph viewer

self.getAvgHandles(); % create measures
s = SolverJMT(self, Solver.defaultOptions, jmtGetPath);
s.writeJSIM(self);
jsimFile = [s.getFilePath, filesep, s.getFileName, '.jsim'];
jsimgView(jsimFile);
% viewerPath = lineViewerGetPath();
% system(sprintf('java -jar "%s" "%s" &', viewerPath, jsimFile));
end
