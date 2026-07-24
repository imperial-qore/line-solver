function jsimwView(self)
% JSIMWVIEW()

[Q,U,R,T,A] = getAvgHandles(self); % create measures
s=SolverJMT(self,Solver.defaultOptions,jmtGetPath);
s.jsimwView;
end