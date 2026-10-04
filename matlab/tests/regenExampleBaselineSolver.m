function regenExampleBaselineSolver(exampleName, solverName)
% REGENEXAMPLEBASELINESOLVER Refresh ONLY one solver's rows of an example baseline.
%
% genAllTestsExamples re-records the whole testsExamples/<example>.mat, so every
% solver in it is re-baselined at once and an unrelated drift is absorbed
% unseen. This refreshes the entries of ONE solver, leaving the rest of the
% file byte-identical, which is what a solver-scoped algorithm change needs.
%
% Usage (from matlab/, after lineStart):
%   (no sibling checkout needed: matlab/tests/ is on the path after lineStart)
%   regenExampleBaselineSolver('cache_replc_fifo', 'SolverSSA')
%
% Every replaced entry is reported with the max abs change per metric, so a
% silent re-baseline is impossible.

cwd = fileparts(mfilename('fullpath'));
baselineFile = fullfile(cwd, 'testsExamples', [exampleName, '.mat']);
if ~isfile(baselineFile)
    line_error(mfilename, sprintf('No baseline at %s', baselineFile));
end
saved = load(baselineFile);
AvgTableEx = saved.AvgTable;

% Re-run the example exactly as the generator does; evalc keeps its own
% per-solver table printouts out of this report.
evalc(exampleName);

fields = {'QLen','Util','RespT','Tput'};
nRepl = 0;
for s = 1:numel(solver)
    if ~strcmp(solver{s}.getName(), solverName)
        continue;
    end
    for f = 1:numel(fields)
        % IndexedTable exposes no Properties, so probe the column the way
        % runTestExample reads it and skip the metrics this table lacks.
        try
            oldv = AvgTableEx{s}.(fields{f})(:);
            newv = AvgTable{s}.(fields{f})(:);
        catch
            continue;
        end
        keep = ~isnan(oldv);
        if numel(oldv) == numel(newv)
            delta = max(abs(newv(keep) - oldv(keep)));
        else
            delta = Inf;
        end
        fprintf('  %-28s %-10s %-6s max|d| = %.3e\n', exampleName, solverName, fields{f}, delta);
    end
    AvgTableEx{s} = AvgTable{s};
    nRepl = nRepl + 1;
end
if nRepl == 0
    line_error(mfilename, sprintf('%s runs no %s', exampleName, solverName));
end

saved.AvgTable = AvgTableEx;
save(baselineFile, '-struct', 'saved');
fprintf('=== %s: %d %s entries refreshed ===\n', exampleName, nRepl, solverName);
end
