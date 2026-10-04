evalc(exampleName);
cwd = fileparts(which(mfilename));

% Under a non-MATLAB backend the baselines describe a different implementation,
% see LINETESTSMOKEMODE: assert that every solver completed and returned a
% usable AvgTable, and leave the recorded values alone.
if lineTestSmokeMode()
    smokeFields = {'QLen','Util','RespT','Tput'};
    smokeNames = cell(1,length(solver));
    for s = 1:length(solver)
        smokeNames{s} = solver{s}.getName;
    end
    smokeValues = cell(1,numel(smokeFields));
    for f = 1:numel(smokeFields)
        vals = cell(1,length(solver));
        for s = 1:length(solver)
            try
                vals{s} = AvgTable{s}.(smokeFields{f})(:);
            catch ME
                error('lineTest:smoke', '%s: %s has no %s column in AvgTable (%s).', ...
                    exampleName, smokeNames{s}, smokeFields{f}, ME.message);
            end
        end
        smokeValues{f} = vals;
    end
    lineTestSmokeCheckSet(exampleName, smokeNames, smokeFields, smokeValues);
    return
end

baselineFile = [cwd,'/testsExamples/',exampleName,'.mat'];
saved = load(baselineFile,'AvgTable');
AvgTableEx = saved.AvgTable;
fields = {'QLen','Util','RespT','Tput'};
for s=1:length(solver)
    isSolverDES = strcmp(solver{s}.getName, 'SolverDES') || strcmp(solver{s}.getName, 'SolverLDES');
    tol = GlobalConstants.FineTol;
    relTol = 0;
    % SolverDES/SolverLDES IS A SEEDED SAMPLE PATH, so its gate has to be the
    % Monte Carlo error the run actually carries, not machine precision.
    % Measured on mqn_singleserver_ps, five seeds at 1e5 samples span 5.8% to
    % 6.8% per row -- the same band as the SolverJMT entry below. An ABSOLUTE
    % 1e-3 on a QLen of 98.57 instead asks for FIVE significant figures of a
    % path whose own seed-to-seed sd is 0.0248, which tests bit-reproduction
    % and not accuracy. That is not hypothetical: routing the C++ engine's
    % variates through fdlibm log1p -- needed because glibc 2.35 and 2.39
    % return different last bits, so ONE binary gave two sample paths --
    % moved Q1/c2 by 0.0013 on a row whose sd across seeds is 0.9386, and the
    % example failed. As with SolverFLD below, the RELATIVE slope moves and
    % the CoarseTol floor stays, so a small station is still held to 1e-3 and
    % an exact zero to 1e-3; tolPrint stays at FineTol, so a drift is PRINTED
    % before it is tolerated. See _kb/08-build-and-test.md.
    if isSolverDES
        tol = GlobalConstants.CoarseTol;
        relTol = 0.10;
    end
    % A SolverJMT row is a seeded sample path recorded at the example's own
    % sample count, usually 1e4. JMT reproduces exactly for a given .jsim, so
    % an edit to the EXPORTER alone moves it while the model stays the same;
    % that happened on 2026-08-14 and moved 15 fork-join rows at once. The
    % assertion gate is therefore relative and sized to the Monte Carlo error
    % those runs actually carry -- five seeds at 1e4 span 4.9% to 8.3% on the
    % fork-join examples. Any movement at all is still PRINTED (relTolPrint
    % below stays at FineTol), so a drift is visible before it is tolerated.
    % See _kb/08-build-and-test.md.
    if strcmp(solver{s}.getName, 'SolverJMT')
        tol = GlobalConstants.CoarseTol;
        relTol = 0.10;
    end
    % SolverLQNS AND SolverLQSIM ARE HELD TO WHAT lqns PRINTS. The wrapper
    % parses the .lqxo, and lqns writes its attributes at SIX SIGNIFICANT
    % DIGITS (service-time="3.22581" on 11-interlock.lqnx), so a parsed value
    % carries a rounding residue of half a quantum, 5e-6 relative. Since
    % 05ef92ac3 an ENTRY's proc-utilization is SUMMED over its activity graph
    % -- lqns reports 0 for an entry in the activity-graph form, which is the
    % only form writeXML emits -- so the residues of a few printed values add.
    % Against the ABSOLUTE FineTol that surfaced as a failure no LINE change
    % can fix: lqn_serial moved one Util row by 5.0e-7 (0.5327135 summed
    % against lqns' own 0.532713), i.e. two orders of magnitude finer than the
    % binary prints. The gate is therefore relative and sized to that
    % printing, not to LINE's arithmetic. A SolverLQNS on a flat Network (its
    % qns* methods, formerly SolverQNS) keeps the absolute gate it always had.
    if (strcmp(solver{s}.getName, 'SolverLQNS') && isa(solver{s}.model, 'LayeredNetwork')) || strcmp(solver{s}.getName, 'SolverLQSIM')
        relTol = 1e-5;
    end
    % SOLVERFLD IS HELD RELATIVELY FOR THE SAME KIND OF REASON. Its answer is
    % the end state of an ODE integration run at AbsTol = RelTol = options.tol,
    % which defaults to 1e-4; an ABSOLUTE 1e-8 on a metric of magnitude 13.57
    % asks for thirteen significant digits from it, and what actually decides
    % the last ones is the host's BLAS. init_state_ps failed on picard04 at
    % 6.6e-9 on RespT = 13.5703845 and passes on picard09, against code that is
    % the same on both. Reading FineTol RELATIVELY keeps the same gate on a
    % metric of order one and stops a large one from demanding accuracy the
    % integrator never claimed; the absolute floor below still applies.
    %
    % 1e-5, NOT FineTol, BECAUSE THE MOMENT CLOSURE CANNOT PROMISE MORE. The
    % second-order methods ('minnormal' and the 'default' that resolves to it)
    % alternate a mean solve against a Lyapunov covariance until sigma2 settles
    % to SOLVER_FLUID_MOMENTS' own mom_tol, 1e-6 relative. What comes out is
    % reproducible to about that, not to 1e-8: measured on mqn_singleserver_ps,
    % picard03 against picard06, same commit and same R2026a, the deviation is
    % a near-CONSTANT 7.19e-7 RELATIVE across every cell it touches --
    %
    %   QLen  Queue2/Closed  0.596208   |d| 4.29e-07   rel 7.19e-07
    %   QLen  Queue2/Open    0.361353   |d| 2.60e-07   rel 7.19e-07
    %   RespT Queue2/Closed  0.851724   |d| 6.13e-07   rel 7.19e-07
    %   RespT Queue2/Open    1.204515   |d| 8.66e-07   rel 7.19e-07
    %
    % -- a uniform offset in the fixed point the two hosts converge to, not
    % scattered noise. It surfaces on the SMALL stations because there
    % max(tol, relTol*|v|) collapses onto the 1e-8 ABSOLUTE floor, which asks
    % eight decimals of a queue length of 0.36 while asking the 98.9 row for
    % only 9.9e-7. So the RELATIVE slope moves and the floor stays: an exact
    % zero is still held to 1e-8. 1e-5 leaves an order of magnitude over the
    % measured 7.2e-7 and still gates five significant digits, which no real
    % modelling regression survives. Same reasoning, and the same value, as the
    % SolverLQNS entry above. tolPrint stays at FineTol, so a drift is still
    % PRINTED before it is tolerated.
    if strcmp(solver{s}.getName, 'SolverFLD') || strcmp(solver{s}.getName, 'SolverFluid')
        relTol = 1e-5;
    end

    for f=1:length(fields)
        try
            val = AvgTable{s}.(fields{f})(:);
            valex = AvgTableEx{s}.(fields{f})(:);
        catch ME
            printSolverDiagnostics(solver{s}, fields{f}, exampleName, tol, baselineFile, s, AvgTableEx{s});
            fprintf('Failed while reading AvgTable entries: %s\n', ME.message);
            rethrow(ME);
        end
        val(isnan(valex)) = [];
        valex(isnan(valex)) = [];
        % tolv is what the assertion uses; tolPrint is what makes a movement
        % visible, and stays tight even where the assertion is deliberately loose
        tolv = max(tol, relTol*abs(valex(:)));
        tolPrint = GlobalConstants.FineTol;
        if any(abs(val(:)-valex(:))>tolPrint)
            printSolverDiagnostics(solver{s}, fields{f}, exampleName, tol, baselineFile, s, AvgTableEx{s});
            if length(val)==length(valex)
                fprintf(sprintf('[%s: New values (%s), Old values, Diff]:\n',solver{s}.name,fields{f}));
                [val,valex,val-valex]
            else
                fprintf(sprintf('%s: New values (%s):\n',solver{s}.name,fields{f}));
                val
                fprintf('Old values:\n');
                valex
            end
        end

        if any(isinf(val))
            fprintf('Warning: %s returned an Inf entry in AvgTable. Skipping test. \n',solver{s}.getName);
        else
            if ~all(abs(val-valex)<tolv)
                printSolverDiagnostics(solver{s}, fields{f}, exampleName, tol, baselineFile, s, AvgTableEx{s});
                if isa(model,'Network')
                    if SolverJMT.supports(model)
                        fprintf('JMT avgTable:\n');
                        SolverJMT(model,'samples',1e6).getAvgTable
                    end
                end
            end
            assert(all(abs(val-valex)<tolv),[solver{s}.getName,' failed on ',fields{f},' in ',exampleName,'.m']);
        end
    end
end

function printSolverDiagnostics(solverEntry, fieldName, exampleName, tol, baselineFile, baselineIndex, baselineEntry)
fprintf('Diagnostic: example=%s field=%s tol=%.12g baseline=%s\n', ...
    exampleName, fieldName, tol, baselineFile);
fprintf('Diagnostic: solver=%s solverClass=%s\n', ...
    getSolverName(solverEntry), class(solverEntry));

global LINEDefaultLang
if isempty(LINEDefaultLang)
    fprintf('Diagnostic: LINEDefaultLang=<empty>\n');
else
    fprintf('Diagnostic: LINEDefaultLang=%s\n', char(string(LINEDefaultLang)));
end

if isprop(solverEntry,'options') && isstruct(solverEntry.options)
    opts = solverEntry.options;
    fprintf('Diagnostic: options.lang=%s options.method=%s options.seed=%s options.samples=%s options.verbose=%s\n', ...
        getStructFieldAsText(opts,'lang'), ...
        getStructFieldAsText(opts,'method'), ...
        getStructFieldAsText(opts,'seed'), ...
        getStructFieldAsText(opts,'samples'), ...
        getStructFieldAsText(opts,'verbose'));
    if isfield(opts,'config') && isstruct(opts.config)
        fprintf('Diagnostic: options.config.eventcache=%s options.config.multiserver=%s options.config.fork_join=%s\n', ...
            getStructFieldAsText(opts.config,'eventcache'), ...
            getStructFieldAsText(opts.config,'multiserver'), ...
            getStructFieldAsText(opts.config,'fork_join'));
    end
end

fprintf('Diagnostic: baseline AvgTableEx{%d} class=%s\n', baselineIndex, class(baselineEntry));
if isstruct(baselineEntry)
    fprintf('Diagnostic: baseline fields=%s\n', strjoin(fieldnames(baselineEntry), ', '));
end
end

function out = getSolverName(solverEntry)
out = class(solverEntry);
try
    out = solverEntry.getName();
catch
    if isprop(solverEntry,'name')
        out = solverEntry.name;
    end
end
out = char(string(out));
end

function out = getStructFieldAsText(s, fieldName)
if ~isfield(s, fieldName)
    out = '<missing>';
    return
end

value = s.(fieldName);
if ischar(value) || isstring(value)
    out = char(string(value));
elseif isnumeric(value) || islogical(value)
    out = mat2str(value);
elseif isempty(value)
    out = '<empty>';
else
    out = class(value);
end
end
