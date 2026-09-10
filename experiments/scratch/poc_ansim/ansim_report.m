function ansim_report(R, csvpath)
% ANSIM_REPORT  Print the comparison table and write it as CSV.
% R is a cell array of per-configuration result structs from poc_ansim_prio.

methods = {'dpsgap','holgap','prsgap','fcfsgap', ...
           'ansim_ctmcP','ansim_live','ansim_mn', ...
           'FLDmn','MVAdps','NCdps', ...
           'MVAcl','MVAshadow','MAMhol', ...
           'MVAprs','MAMprs','MVAmarie', ...
           'LDES10000','LDES100000','SSA'};

fid = fopen(csvpath, 'w');
fprintf(fid, 'law,scv,N');
for j = 1:numel(methods), fprintf(fid, ',err_%s', methods{j}); end
for j = 1:numel(methods), fprintf(fid, ',t_%s', methods{j}); end
fprintf(fid, '\n');

fprintf('\n\n%-10s %5s %3s', 'service', 'SCV', 'N');
for j = 1:numel(methods), fprintf(' %11s', methods{j}); end
fprintf('\n');

for i = 1:numel(R)
    r = R{i};
    fprintf('%-10s %5.1f %3d', r.law, r.scv, r.N);
    fprintf(fid, '%s,%.1f,%d', r.law, r.scv, r.N);
    for j = 1:numel(methods)
        v = getf(r, ['err_' methods{j}]);
        if isnan(v), fprintf(' %11s', '--'); else, fprintf(' %11.4f', v); end
        fprintf(fid, ',%g', v);
    end
    for j = 1:numel(methods)
        fprintf(fid, ',%g', getf(r, ['t_' methods{j}]));
    end
    fprintf('\n'); fprintf(fid, '\n');
end

% Timing table, geometric-mean seconds per method over the configurations that
% produced a number.
fprintf('\n%-12s %10s %10s\n', 'method', 'mean L1', 'mean sec');
for j = 1:numel(methods)
    e = arrayfun(@(i) getf(R{i}, ['err_' methods{j}]), 1:numel(R));
    t = arrayfun(@(i) getf(R{i}, ['t_'   methods{j}]), 1:numel(R));
    e = e(~isnan(e)); t = t(~isnan(t));
    if isempty(e)
        fprintf('%-12s %10s %10s\n', methods{j}, '--', '--');
    else
        fprintf('%-12s %10.4f %10.2f\n', methods{j}, mean(e), mean(t));
    end
end

fclose(fid);
fprintf('\nwrote %s\n', csvpath);
end

function v = getf(s, f)
if isfield(s, f) && isnumeric(s.(f)) && isscalar(s.(f))
    v = s.(f);
else
    v = NaN;
end
end
