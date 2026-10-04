testName = 'tut06_cache_lru_zipf';
evalc(testName);
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end

TOL=1e-2; % 1-percent tolerance

% The default SSA back-end is the NRM engine in every codebase, but the
% random streams are not shared: MATLAB draws from its own rand/randn, while
% the JAR (commons-math3 MersenneTwister via RandomManager) and native Python
% draw from an identical MT19937 stream and agree to ~1e-8. At seed 23000 and
% samples 2e4 the two streams are therefore distinct, equally valid
% realizations (hit probability 0.8313 vs 0.8241), and they stay ~0.8% apart
% at samples 5e5, i.e. above TOL. The baseline is consequently keyed on
% options.lang. Each branch is exactly reproducible (verified over 3 repeats
% per language).
global LINEDefaultLang %#ok<GVMIS>
if ~isempty(LINEDefaultLang) && any(strcmp(LINEDefaultLang, {'java','python'}))
    Qsaved = [
                       0
       0.483812287110524
       0.516187594651938];
    Tsaved = [
       2.99565777672989
       2.41906143555261
       0.516187594651938];
else
    % Re-recorded 2026-08-25 against the corrected NRM engine (line-dev
    % a07778949: BLOCKED departures are subtracted from the propensity
    % integral, and utilization is work-based). That changed the MATLAB sample
    % path deterministically -- three runs in one session agree to 15 digits --
    % and moved the hit probability from 0.8313 to 0.8260, i.e. TOWARD the
    % java/python branch above (0.8241), which is what running the same engine
    % should do. The previous values were recorded 2026-07-22, before it.
    Qsaved = [
                       0
       0.487003306430337
       0.512996576087011];
    Tsaved = [
       2.94881037925967
       2.43501653215169
       0.512996576087011];
end
Usaved = Qsaved; % Util == QLen: both delays are infinite-server stations

QN = ssaAvgTable.QLen(:);
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = ssaAvgTable.Util(:);
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = ssaAvgTable.RespT(:);
Rsaved = [
                   0
   0.200000000000000
   1.000000000000000];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ssaAvgTable.Tput(:);
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

lineTestAssertChecks(testName, gsFailures);
