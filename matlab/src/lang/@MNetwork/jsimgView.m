function jsimgView(self)
% JSIMGVIEW()

self.getAvgHandles(); % create measures
s=SolverJMT(self,Solver.defaultOptions,jmtGetPath); s.jsimgView;
end