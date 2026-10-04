function results = lineTestAccum(results)
% LINETESTACCUM Record runtests results so a failing suite fails the process.
%
%   RESULTS = LINETESTACCUM(RESULTS) appends RESULTS to the suite-wide
%   accumulator and passes them through unchanged.
%
%   runtests returns a TestResult array rather than raising, so a script that
%   discards the return value always leaves MATLAB with exit status 0: a suite
%   with errored tests then reports success to the caller (the "All tests
%   completed!" false green). Every allTests*.m routes its runtests output
%   through here, and lineTestAssert raises at the end of allTests, which makes
%   matlab -batch exit nonzero.
%
%   HELD IN ROOT APPDATA, NOT A GLOBAL (2026-09-09). Thirteen examples in the
%   corpus run `clear all`, which empties the global workspace -- so an
%   accumulator kept in a global loses every result recorded before the last
%   such example. That is not hypothetical: run 20260908_230513 reported
%   `Suite totals: 209 tests, 0 failed` and `All tests passed.` with exit code
%   0, where 209 is the size of allTestsParity ALONE, the final block. The
%   thousands of results before it -- including the errored
%   allTestsCSRunner/test_CQN_Cox_CS_1, whose ten JMT mismatches are printed in
%   the same log -- had been wiped, so the phase reported success over a red
%   test. Root appdata survives `clear all`, which is why ParityRows.globals()
%   already keeps the LINE global snapshot there for the same reason.
%
%   See also LINETESTASSERT, LINETESTRESET.

prev = getappdata(0, 'LINETestResults');
if isempty(prev)
    setappdata(0, 'LINETestResults', results);
else
    setappdata(0, 'LINETestResults', [prev, results]);
end
end
