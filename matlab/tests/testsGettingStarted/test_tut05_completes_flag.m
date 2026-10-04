evalc('tut05_completes_flag');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end

TOL=1e-2; % 1-percent tolerance

QN = ncAvgTable.QLen(:);
Qsaved = [
    0.1667
    0.3333
    0.5000];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = ncAvgTable.Util(:);
Usaved = [
    0.1667
    0.3333
    0.5000];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = ncAvgTable.RespT(:);
Rsaved = [
    1
    2
    3];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ncAvgTable.Tput(:);
Tsaved = [
    0.1667
    0.1667
    0.1667];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

RN = ncAvgSysTable.SysRespT(:);
Rsaved = [2];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ncAvgSysTable.SysTput(:);
Tsaved = [0.5];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

RN = ncAvgSysTable2.SysRespT(:);
Rsaved = [6];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ncAvgSysTable2.SysTput(:);
Tsaved = [0.16667];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

lineTestAssertChecks(testName, gsFailures);
