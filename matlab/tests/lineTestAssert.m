function lineTestAssert()
% LINETESTASSERT Raise if any accumulated test failed or errored.
%
%   Turns the accumulated runtests results into a process exit status: under
%   matlab -batch an uncaught error yields a nonzero exit, which is what lets
%   run_tests.sh report the suite honestly.
%
%   See also LINETESTACCUM, LINETESTRESET.

% Root appdata rather than a global: a global does not survive the `clear all`
% that thirteen corpus examples run, and losing the accumulator here is silent
% -- it reports a short, all-green suite. See LINETESTACCUM.
LINETestResults = getappdata(0, 'LINETestResults');

if isempty(LINETestResults)
    fprintf('lineTestAssert: no test results were recorded.\n');
    return
end

failed = [LINETestResults.Failed];
incomplete = [LINETestResults.Incomplete];
bad = failed | incomplete;

fprintf('\n==========================================\n');
fprintf('Suite totals: %d tests, %d failed, %d incomplete\n', ...
    numel(LINETestResults), sum(failed), sum(incomplete));

if ~any(bad)
    fprintf('All tests passed.\n');
    fprintf('==========================================\n');
    return
end

names = {LINETestResults(bad).Name};
fprintf('Failing tests:\n');
for i = 1:numel(names)
    fprintf('  %s\n', names{i});
end
fprintf('==========================================\n');

error('lineTestAssert:suiteFailed', '%d of %d tests failed or did not complete.', ...
    sum(bad), numel(LINETestResults));
end
