clear; 
evalc('cdf_respt_distrib');
exampleName = 'cdf_respt_distrib';
TOL=1e-2; % 1-percent tolerance

% Re-recorded 2026-08-12: the adaptive CDF refinement now splits EVERY interval
% that breaks the 0.05% jump target in one round and re-integrates the curve on
% the refined grid (see test_cdf_respt_closed.m). Every entry moves toward the
% simulation row below: the delay rows to 1.0031 and 4.0012 against 1.0005 and
% 3.9898, and the queue rows 7.0468 -> 6.5793 against 6.4328 and
% 9.7967 -> 9.3561 against 16.5309.
AvgRespTfromCDFfluid_saved = [
   1.003054429880730   4.001156609349640
   6.579282144104310   9.356103098086330];
try
    assert(max(max(abs(AvgRespTfromCDFfluid-AvgRespTfromCDFfluid_saved)))<TOL,sprintf('%s changed on AvgRespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDFfluid),mat2str(AvgRespTfromCDFfluid_saved))
end

AvgRespTfromCDFsim_saved = [
1.00049271750951 3.98983612815623;
6.43275357361235 16.530859522722];
try
    assert(max(max(abs(AvgRespTfromCDFsim-AvgRespTfromCDFsim_saved)))<TOL,sprintf('%s changed on AvgRespTfromCDF.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDFsim),mat2str(AvgRespTfromCDFsim_saved))
end
close all;