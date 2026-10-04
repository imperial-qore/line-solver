clear; 
evalc('cdf_respt_closed_threeclasses');
exampleName = 'cdf_respt_closed_threeclasses';
TOL=1e-2; % 1-percent tolerance

AvgRespT_saved = [
     1     1     0
     1     4     0];
try
    assert(max(max(abs(AvgRespT-AvgRespT_saved)))<TOL,sprintf('%s changed on AvgRespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespT),mat2str(AvgRespT_saved))
end

% Re-recorded 2026-08-12: the adaptive CDF refinement now splits EVERY interval
% that breaks the 0.05% jump target in one round and re-integrates the curve on
% the refined grid. This model used to come back on 130 points over [0,200] and
% its right-endpoint mean read 1.1265 for a response time that is exactly
% Exp(1); it now reads 1.0029 against the AvgRespT row above.
AvgRespTfromCDF_saved = [
   1.002874300199440   1.002874300199440                   0
   1.002874300199440   4.003443315167100                   0];
try
    assert(max(max(abs(AvgRespTfromCDF-AvgRespTfromCDF_saved)))<TOL,sprintf('%s changed on AvgRespTfromCDF.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDF),mat2str(AvgRespTfromCDF_saved))
end

% Re-recorded 2026-08-12 with the same refinement: the exponential rows land on
% 1.0101 against the exact 1.0, where the coarse grid read 0.8690.
SqCoeffOfVariationRespTfromCDF_saved = [
   1.010062517359530   1.010062517359530                 NaN
   1.010062517359530   0.503899386442527                 NaN];
try
    assert(max(max(abs(SqCoeffOfVariationRespTfromCDF-SqCoeffOfVariationRespTfromCDF_saved)))<TOL,sprintf('%s changed on SqCoeffOfVariationRespTfromCDF.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(SqCoeffOfVariationRespTfromCDF),mat2str(SqCoeffOfVariationRespTfromCDF_saved))
end