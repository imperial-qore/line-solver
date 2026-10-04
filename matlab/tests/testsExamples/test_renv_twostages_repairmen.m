clear; 
exampleName = 'renv_twostages_repairmen';
evalc('renv_twostages_repairmen');
TOL=1e-2;
QN = AvgTable.QLen(:);
Qsaved = [
    0.5591
    0.4409];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTable.Util(:);
Usaved = [
    0.5591
    0.4409];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

TN = AvgTable.Tput(:);
Tsaved = [
    0.6946
    0.6838];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
