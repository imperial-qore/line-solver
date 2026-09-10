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

if strcmp(options.lang,'cpp')
    line_debug(options, '%s: using lang=cpp, delegating to the C++ line-cli', tag);
    Solver.resetRandomGeneratorSeed(options.seed);
    [QN,UN,RN,TN,AN,WN,runtime] = CPPLINE.getAvg(self.name, self.model, options);
    self.setAvgResults(QN,UN,RN,TN,AN,WN,[],[],runtime,options.method,NaN);
    handled = true;
    return
end

Solver.resetRandomGeneratorSeed(options.seed);

% MODEL TRANSFORMATION, opt-in through options.config.transform. The strategy
% rewrites the model into subproblems, TRANSFORMSOLVE solves each with an
% instance of THIS solver and maps the metrics back, so a transformation
% written once serves every solver that reaches this preamble rather than the
% one it was first written for.
%
% This is the seam because it is the one opening every runAnalyzer already
% shares, and it already carries a HANDLED early-return contract. Wiring the
% dispatch here reaches all eight of its callers (MVA, NC, CTMC, FLD, SSA, MAM,
% AG, BA) without editing any of them. The JMT, QNS and LDES wrappers inherit
% transformSolve but do not call this preamble, so a transform method name is ignored
% there rather than honoured.
%
% The inner solve carries transform='none', so a transformed submodel arriving
% back here cannot re-enter the driver.
if isfield(options,'config') && isstruct(options.config) ...
        && isfield(options.config,'transform') && ~isempty(options.config.transform) ...
        && ~strcmpi(options.config.transform,'none')
    line_debug(options, '%s: model transformation ''%s''', tag, options.config.transform);
    tr = transformSolve(self, options);
    runtime = tr.runtime;
    sn = getStruct(self);
    T = getAvgTputHandles(self);
    AN = sn_get_arvr_from_tput(sn, tr.TN, T);
    self.setAvgResults(tr.QN, tr.UN, tr.RN, tr.TN, AN, [], tr.CN, tr.XN, ...
        runtime, [options.method '/' tr.method], tr.iter);
    handled = true;
    return
end
end
