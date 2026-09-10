function lineClearWorkspace(maxAgeHours)
% LINECLEARWORKSPACE Remove stale temporary files from past LINE executions.
%
% The workspace folder (tempdir/line_workspace) is a process-global path
% shared by every LINE process on the machine. Solvers create per-run temp
% subdirectories there (e.g. JMT writes model.jsim / .jsim-result.jsim, CTMC
% saves the generator .mat). Deleting a subdirectory that another concurrently
% running MATLAB/JVM is still using makes that solver fail mid-run -- e.g. JMT
% aborts with "java.io.IOException: No such file or directory" in copyFile when
% its result directory disappears.
%
% To stay safe under concurrency (mirroring the JAR, which uses unique temp
% dirs and never bulk-clears the workspace) this only removes entries that have
% not been modified recently: recently-touched entries may belong to a solve
% running in another process. Genuinely stale leftovers from crashed prior runs
% are older than the threshold and still get cleaned. MAXAGEHOURS defaults to 1,
% far above any realistic single-solve duration in the example/parity suites.

if nargin < 1 || isempty(maxAgeHours)
    maxAgeHours = 1;
end
staleAgeDays = maxAgeHours / 24;

workspaceFolder = fullfile(tempdir, 'line_workspace/');
if ~exist(workspaceFolder, 'dir')
    return;
end
files = dir(workspaceFolder);
nowNum = now;
for k = 1:length(files)
    name = files(k).name;
    if ismember(name, {'.', '..'})
        continue;
    end
    % Skip entries touched within the stale window: they may be in use by a
    % LINE solve running concurrently in another process.
    if (nowNum - files(k).datenum) < staleAgeDays
        continue;
    end
    target = fullfile(workspaceFolder, name);
    if isfolder(target)
        rmdir(target, 's');
    else
        delete(target);
    end
end
end
