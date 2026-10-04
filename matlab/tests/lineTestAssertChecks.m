function lineTestAssertChecks(testName, failures)
% LINETESTASSERTCHECKS Fail the calling test if any of its soft checks failed.
%
%   LINETESTASSERTCHECKS(TESTNAME, FAILURES) raises 'lineTest:checksFailed'
%   when the cell array FAILURES is non-empty, and returns silently otherwise.
%
%   A getting-started test runs several baseline checks, each wrapped in
%   try/catch so that every mismatch is PRINTED rather than only the first.
%   Until 2026-10-04 the catch was the end of it: the assertion was swallowed,
%   the test passed, and runtests counted it green. Run 20261003_210029 printed
%   "Assertion failed" for tut01, tut02, tut04 and tut08 while its Failure
%   Summary and `Suite totals: 1045 tests, 4 failed` named only examples. So
%   each catch now appends me.message to FAILURES, and the test ends by calling
%   this, which turns the collected mismatches into a real test failure.
%
%   Same shape as runTest.m's checkAssert/error pair in line-test.git.
%
%   See also LINETESTACCUM, LINETESTASSERT.

if isempty(failures)
    return
end
error('lineTest:checksFailed', '%s: %d check(s) failed: %s', ...
    testName, numel(failures), strjoin(failures, ' | '));
end
