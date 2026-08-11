function [runtime, sruntime, results] = iterate(self, options)
% [RUNTIME, SRUNTIME, RESULTS] = ITERATE()
T0 = tic;
it = 0;
options = self.options;
E = getNumberOfModels(self);
results = cell(1,E);
sruntime = zeros(1,E); % solver runtimes
init(self);
% nearly identical, but parfor based
while ~self.converged(it) && it < options.iter_max
    it = it + 1;
    line_debug('EnsembleSolver iteration %d starting (max=%d)', it, options.iter_max);
    self.pre(it);
    sruntime(it,1:E) = 0;
    T1=tic;
    switch options.method
        case {'parallel'}
            parfor e = self.list(it)
                [results{it,e}, solverTime] = self.analyze(it,e);
                sruntime(it,e) = sruntime(it,e) + solverTime;
            end
        otherwise
            for e = self.list(it)
                line_debug('Analyzing ensemble model %d at iteration %d', e, it);
                [results{it,e}, solverTime] = self.analyze(it,e);
                sruntime(it,e) = sruntime(it,e) + solverTime;
            end
    end
    self.results = results;
    if options.verbose
        Tsolve(it)=toc(T1);
        Ttot=toc(T0);
        line_printf('\nIter %2d. ',it);
    end
    T2=tic;
    self.post(it);
    Tsynch(it)=toc(T2);
    if options.verbose
        line_printf('\b Analyze time: %.3fs. Update time: %.3fs. Runtime: %.3fs. ',Tsolve(it),Tsynch(it),Ttot);
    end
end
finish(self);
runtime = toc(T0);
if options.verbose
    line_printf('\nSummary: Analyze avg time: %.3fs. Update avg time: %.3fs. Total runtime: %.3fs. ',mean(Tsolve),mean(Tsynch),runtime);
end
end