function tf = lineTestSmokeMode()
% LINETESTSMOKEMODE  True when the suite must assert COMPLETION, not goldens.
%
% The recorded values (test/regression/*.mat, test/testsExamples/*.mat) say what
% the MATLAB implementation computes. A DISPATCH run -- lang='java', 'python' or
% 'cpp' -- is a different implementation answering the same model, so its numbers
% legitimately differ, and comparing them against those files reports an engine
% difference as a regression. Under such a backend runTest.m and runTestExample.m
% switch to SMOKE MODE instead: every solver must complete and return usable
% numbers, and the recorded values are neither read nor overwritten.
%
% The backend is taken from the LINEDefaultLang global that the wrapper runners
% set, with the LINE_DEFAULT_LANG environment variable as the durable fallback: a
% test script running `clear all` DELETES the global while the environment
% survives, which is also why SolverOptions.m re-seeds one from the other.
%
% LINE_TEST_SMOKE overrides both directions -- '0' restores the golden comparison
% under a dispatch backend, '1' asks for smoke mode under MATLAB itself.
%
% See matlab/run_tests_wrapper_java.sh and _kb/08-build-and-test.md.

override = strtrim(getenv('LINE_TEST_SMOKE'));
if ~isempty(override)
    tf = ~any(strcmpi(override, {'0', 'false', 'no', 'off'}));
    return
end

global LINEDefaultLang

lang = '';
if ~isempty(LINEDefaultLang) && (ischar(LINEDefaultLang) || isstring(LINEDefaultLang))
    lang = char(LINEDefaultLang);
end
if isempty(lang)
    lang = getenv('LINE_DEFAULT_LANG');
end

tf = ~isempty(lang) && ~strcmpi(strtrim(lang), 'matlab');
end
