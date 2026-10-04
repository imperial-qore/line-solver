clear; 
evalc('cdf_respt_open_twoclasses');
exampleName = 'cdf_respt_open_twoclasses';
TOL=1e-2; % 1-percent tolerance

AvgRespTfromCDFfluid_saved = [
                   0                   0
   1.097874431999490   1.097874431999490
   1.078350188965910   1.078350189590920
   ];
try
    assert(max(max(abs(AvgRespTfromCDFfluid-AvgRespTfromCDFfluid_saved)))<TOL,sprintf('%s changed on AvgRespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDFfluid),mat2str(AvgRespTfromCDFfluid_saved))
end

AvgRespTfromCDFsim_saved = [
                   0                   0
   1.935842575118770   1.950174785401700
   2.090394050286220   2.068682104526130
   ];
try
    assert(max(max(abs(AvgRespTfromCDFsim-AvgRespTfromCDFsim_saved)))<TOL,sprintf('%s changed on AvgRespTfromCDF.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDFsim),mat2str(AvgRespTfromCDFsim_saved))
end
