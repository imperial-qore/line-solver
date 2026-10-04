evalc('tut01_mm1_basics');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end
TOL=1e-2; % 1-percent tolerance

QN = AvgTable.QLen(:);
Qsaved = [
         0
    1.0156];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTable.Util(:);
Usaved = [
         0
    0.5066];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgTable.RespT(:);
Rsaved = [
         0
    1.0159];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = AvgTable.Tput(:);
Tsaved = [
    1.0000
    1.0059];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

lineTestAssertChecks(testName, gsFailures);
