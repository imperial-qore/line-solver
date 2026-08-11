function [handled, options, runtime] = runAnalyzerPreamble(self, options, tag)
% [HANDLED, OPTIONS, RUNTIME] = RUNANALYZERPREAMBLE(OPTIONS, TAG)
%
% Common opening of every runAnalyzer: stamp the wall-clock budget marker,
% delegate to native python when options.lang says so and seed the RNG.
%
% The library attribution is NOT printed here: it is pull-based, like Sage's
% sage.misc.citation.get_systems. Call solver.getLibrariesUsed() to obtain the
% list, or solver.showLibraryAttribution() to print it; lineStart reminds the
% user that it exists. The wrapper acknowledgement (line_ack) follows the same
% rule: it prints only at VerboseLevel.DEBUG.
%
% HANDLED is true when the call was served by the python delegation, in which
% case the results are already set and the caller must return RUNTIME
% immediately. TAG is the short solver name used in the debug line (e.g. 'MVA').
%
% Every solver carried its own copy of this block: six copies of the
% lang='python' guard and seven of the attribution print. The RNG seed
% inside the python branch is load-bearing and must stay BEFORE the early
% return: MATLAB-side samplers invoked after the delegated call (sampleSysAggr,
% sample, the getCdfRespT Monte Carlo) draw from MATLAB's stream, so skipping
% it made a seeded example reproducible under lang='matlab' but not under
% lang='python'.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

handled = false;
runtime = 0;
if nargin < 3 || isempty(tag)
    tag = class(self);
end

% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
if ~isfield(options,'timeout_tic') || isempty(options.timeout_tic)
    options.timeout_tic = tic;
end

if strcmp(options.lang,'python')
    line_debug(options, '%s: using lang=python, delegating to native line_solver', tag);
    Solver.resetRandomGeneratorSeed(options.seed);
    [QN,UN,RN,TN,AN,WN,runtime] = PYLINE.getAvg(self.name, self.model, options);
    self.setAvgResults(QN,UN,RN,TN,AN,WN,[],[],runtime,options.method,NaN);
    handled = true;
    return
end

Solver.resetRandomGeneratorSeed(options.seed);
end
