clear; 
exampleName = 'renv_fourstages_repairmen';
evalc('renv_fourstages_repairmen');
TOL=1e-2;
QN = AvgTable.QLen(:);
Qsaved = [
    0.4446
   29.5554];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTable.Util(:);
Usaved = [
    0.4446
    1.0000];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

TN = AvgTable.Tput(:);
Tsaved = [
    0.9780
    1.0000];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
