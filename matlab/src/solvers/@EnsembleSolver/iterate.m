function [runtime, sruntime, results] = iterate(self, options)
% [RUNTIME, SRUNTIME, RESULTS] = ITERATE()
T0 = tic;
it = 0;
options = self.options;
E = getNumberOfModels(self);
results = cell(1,E);
sruntime = zeros(1,E); % solver runtimes
% The status row below is closed on the normal path, after the loop. This
% guard covers the OTHER path: a layer that throws would otherwise leave the
% row open with no newline, and the error message would be backspaced over by
% nothing -- it is simply glued onto the last iteration. CLOSE is idempotent,
% so the normal path is unaffected.
statusGuard = onCleanup(@() LineStatus.close()); %#ok<NASGU>
init(self);
% nearly identical, but parfor based
% The convergence test is hoisted out of the while condition so that the reason
% the loop ended is still known below: exiting on iter_max is a NON-convergence
% and must be reported, not returned as if it were a fixed point. The call
% sequence is unchanged -- converged(0), converged(1), ... converged(n) -- since
% && evaluated it first on every pass anyway.
hasConverged = self.converged(it);
while ~hasConverged && it < options.iter_max
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
    if LineConsole.isActive()
        Tsolve(it)=toc(T1);
        Ttot=toc(T0);
    elseif options.verbose
        % The analyze timings are taken here, but NOTHING IS PRINTED YET -- the
        % row is drawn once, below, when every field on it is a measured value.
        % An in-progress row was tried both ways and neither is worth it: one
        % that shrank to the bare counter blanked a hundred-odd columns for the
        % whole update step, and one that reserved the missing field showed
        % dashes where a number goes, which reads as a measurement that failed
        % rather than one not yet taken. Leaving the PREVIOUS iteration's row
        % standing until this one has real numbers says the same thing and
        % never puts anything on screen that is not a result.
        Tsolve(it)=toc(T1);
        Ttot=toc(T0);
    end
    T2=tic;
    self.post(it);
    Tsynch(it)=toc(T2);
    if LineConsole.isActive()
        Tsolve(it)=toc(T1);
        Ttot=toc(T0);
    elseif options.verbose
        % ONE ROW FOR THE WHOLE RUN, rewritten each iteration rather than one
        % line per iteration: a few hundred iterations otherwise scroll
        % everything else out of the terminal to say the same numbers again.
        % Laid down here and extended by the convergence test, which is what
        % knows the iteration error -- see LineStatus.
        LineStatus.set('Iter %2d. Analyze time: %.3fs. Update time: %.3fs. Runtime: %.3fs.',it,Tsolve(it),Tsynch(it),Ttot);
    end
    % Last thing in the body: converged() extends the status row above with the
    % iteration error, and LineStatus.set would overwrite that if it ran after.
    hasConverged = self.converged(it);
end
if ~hasConverged
    line_warning(mfilename, ['The %s ensemble fixed point did not converge in options.iter_max=%d ' ...
        'iterations; the returned solution is the last iterate and may be far from the fixed point. ' ...
        'Raise options.iter_max, or loosen options.iter_tol only if the residual is already small.\n'], ...
        class(self), options.iter_max);
end
finish(self);
runtime = toc(T0);
if LineConsole.isActive()
    LineConsole.step(['ensemble solved: %.3f s per iteration analyzing, ' ...
        '%.3f s updating, %.3f s in total'], mean(Tsolve), mean(Tsynch), runtime);
elseif options.verbose
    % close the row before the summary: it carries the newline the rewritten
    % row never had, so the summary cannot be glued onto the last iteration.
    LineStatus.close();
    line_printf('Summary: Analyze avg time: %.3fs. Update avg time: %.3fs. Total runtime: %.3fs.\n',mean(Tsolve),mean(Tsynch),runtime);
end
end