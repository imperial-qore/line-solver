function [runtime, analyzer] = runAnalyzer(self, options)
% RUNANALYZER Execute the bound-analysis solver.
%
% Dispatches noniterative bound families to the shared bound handler
% solver_ba_analyzer, and hierarchical/iterative families to their
% SolverBA-native analyzers. Each method returns a single (upper or lower)
% bound; use getBounds() to obtain the {lower,upper} bracket for a family.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = self.getOptions;
end

sn = self.getStruct;

if sn.nclosedjobs <= 0
    line_error(mfilename, 'SolverBA supports closed queueing networks only.');
end

% 'default' -> geometric upper bound; 'qr' -> QRF quadratic reduction; 'lr' ->
% LP linear reduction (bare 'lr' means 'lr.upper'), distinct from
% 'qrf.mmi.linear'. see _kb/06-solver-catalog.md for rationale
method = options.method;
if strcmp(method,'default')
    method = 'gb.upper';
elseif strcmp(method,'lr')
    method = 'lr.upper';
elseif strcmp(method,'qr')
    method = 'qrf.mmi';
end

% Show library attribution (QRF uses the Optimization Toolbox) once.
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    libs = SolverBA.getLibrariesUsed([], setfield(options,'method',method)); %#ok<SFLD>
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

iter = 1;
if startsWith(method,'qrf')
    % QRF (Quadratic Reduction Framework) LP-based bounds for single-class
    % closed networks with PH service. The adapter enforces its own gating.
    analyzer = @(qn) solver_ba_qrf_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,runtime] = solver_ba_qrf_analyzer(sn, bopts);
elseif any(strcmp(method, self.listValidMethods()))
    analyzer = @(qn) solver_ba_analyzer(qn, setfield(options,'method',method)); %#ok<SFLD>
    bopts = options; bopts.method = method;
    [QN,UN,RN,TN,CN,XN,lG,runtime,iter] = solver_ba_analyzer(sn, bopts);
else
    line_error(mfilename, ['Unknown bound method ''%s''. Valid methods: %s'], ...
        method, strjoin(self.listValidMethods(), ', '));
end

M = sn.nstations;
R = sn.nclasses;
% Arrival rates from throughputs via sn_get_arvr_from_tput (not a zero matrix).
% see _kb/06-solver-catalog.md for rationale
AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
WN = zeros(M,R);
self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,method,iter);
end
