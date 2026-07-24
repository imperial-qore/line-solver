function run_all_mexify()
% RUN_ALL_MEXIFY  Run all mexify scripts to generate MEX files.
%
%   RUN_ALL_MEXIFY() runs every mexify.m script found under src/api/ and
%   util/, generating MEX files for performance-critical functions.
%   Each codegen call is wrapped in try-catch so individual failures
%   do not block subsequent compilations.
%
%   Requires MATLAB Coder.

mexify_dirs = { ...
    fullfile('src','api','aoi'), ...
    fullfile('src','api','cache'), ...
    fullfile('src','api','fj'), ...
    fullfile('src','api','mam'), ...
    fullfile('src','api','mc'), ...
    fullfile('src','api','pfqn'), ...
    fullfile('src','api','qsys'), ...
    fullfile('src','api','sn'), ...
    fullfile('util') ...
};

baseDir = fileparts(mfilename('fullpath'));
totalPass = 0;
totalFail = 0;
totalSkip = 0;
failedFunctions = {};

for d = 1:length(mexify_dirs)
    mexDir = fullfile(baseDir, mexify_dirs{d});
    mexFile = fullfile(mexDir, 'mexify.m');
    if ~isfile(mexFile)
        fprintf('[SKIP] %s (not found)\n', mexify_dirs{d});
        continue;
    end
    fprintf('[START] %s\n', mexify_dirs{d});

    % Read mexify.m
    fid = fopen(mexFile, 'r');
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    % Count skipped (commented-out codegen calls)
    dirSkip = 0;
    for li = 1:length(lines)
        line = lines{li};
        if ischar(line) && ~isempty(strtrim(line))
            trimmed = strtrim(line);
            if trimmed(1) == '%' && contains(trimmed, 'codegen') && contains(trimmed, '-config')
                dirSkip = dirSkip + 1;
            end
        end
    end

    % Transform: wrap each codegen call in try-catch, write to temp script
    tmpFile = fullfile(mexDir, 'mexify_wrapped_tmp__.m');
    fid = fopen(tmpFile, 'w');
    % Header: counters
    fprintf(fid, 'dirPass__ = 0; dirFail__ = 0; failedFuncs__ = {};\n');

    blockDepth = 0;
    for li = 1:length(lines)
        line = lines{li};
        if ~ischar(line)
            continue;
        end
        trimmed = strtrim(line);

        % Track nested block comments (%{ ... %})
        if strcmp(trimmed, '%{')
            blockDepth = blockDepth + 1;
            fprintf(fid, '%s\n', line);
            continue;
        end
        if strcmp(trimmed, '%}')
            blockDepth = max(0, blockDepth - 1);
            fprintf(fid, '%s\n', line);
            continue;
        end
        if blockDepth > 0
            fprintf(fid, '%s\n', line);
            continue;
        end

        if isempty(trimmed)
            fprintf(fid, '\n');
            continue;
        end

        % Check if this is an active codegen command
        if startsWith(trimmed, 'codegen ')
            % Extract function name for reporting
            tokens = regexp(trimmed, 'codegen\s+.*?cfg\s+(\w+)', 'tokens');
            if ~isempty(tokens)
                funcName = tokens{1}{1};
            else
                funcName = 'unknown';
            end
            % Wrap in try-catch
            fprintf(fid, 'try\n');
            fprintf(fid, '    %s\n', trimmed);
            fprintf(fid, '    dirPass__ = dirPass__ + 1;\n');
            fprintf(fid, 'catch me__\n');
            fprintf(fid, '    fprintf(''  [FAIL] %s: %%s\\n'', me__.message);\n', funcName);
            fprintf(fid, '    dirFail__ = dirFail__ + 1;\n');
            fprintf(fid, '    failedFuncs__{end+1} = ''%s'';\n', funcName);
            fprintf(fid, 'end\n');
        else
            % Pass through all other lines (setup, comments, etc.)
            fprintf(fid, '%s\n', line);
        end
    end
    fclose(fid);

    % Run the wrapped script
    origDir = pwd;
    cd(mexDir);
    dirPass = 0;
    dirFail = 0;
    dirFailedFuncs = {};
    try
        run('mexify_wrapped_tmp__');
        dirPass = dirPass__;
        dirFail = dirFail__;
        dirFailedFuncs = failedFuncs__;
    catch me
        fprintf('  [ERROR] Script-level error: %s\n', me.message);
    end
    cd(origDir);

    % Clean up temp file
    if isfile(tmpFile)
        delete(tmpFile);
    end
    % Clean up codegen dir if created
    codegenDir = fullfile(mexDir, 'codegen');
    if isfolder(codegenDir)
        rmdir(codegenDir, 's');
    end

    % Report
    if dirFail == 0
        fprintf('[PASS]  %s (%d compiled', mexify_dirs{d}, dirPass);
    else
        fprintf('[DONE]  %s (%d compiled, %d failed', mexify_dirs{d}, dirPass, dirFail);
    end
    if dirSkip > 0
        fprintf(', %d skipped)\n', dirSkip);
    else
        fprintf(')\n');
    end

    for i = 1:length(dirFailedFuncs)
        failedFunctions{end+1} = sprintf('%s/%s', mexify_dirs{d}, dirFailedFuncs{i}); %#ok<AGROW>
    end
    totalPass = totalPass + dirPass;
    totalFail = totalFail + dirFail;
    totalSkip = totalSkip + dirSkip;
end

fprintf('\n=== MEXify Summary ===\n');
fprintf('Compiled: %d\n', totalPass);
if totalFail > 0
    fprintf('Failed:   %d\n', totalFail);
    for i = 1:length(failedFunctions)
        fprintf('  - %s\n', failedFunctions{i});
    end
end
if totalSkip > 0
    fprintf('Skipped:  %d\n', totalSkip);
end
end
