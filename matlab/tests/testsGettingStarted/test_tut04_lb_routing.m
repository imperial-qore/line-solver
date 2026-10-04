evalc('tut04_lb_routing');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end
TOL=1e-2; % 1-percent tolerance

QN = ldesAvgTable.QLen(:);
Qsaved = [
         0
    0.3325
    0.3322];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = ldesAvgTable.Util(:);
Usaved = [
         0
    0.2512
    0.2475];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = ldesAvgTable.RespT(:);
Rsaved = [
         0
    0.6607
    0.6649];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ldesAvgTable.Tput(:);
Tsaved = [
    1.0000
    0.5036
    0.5026];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

QN = ldesAvgTableRR.QLen(:);
Qsaved = [
         0
    0.2883
    0.2840];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = ldesAvgTableRR.Util(:);
Usaved = [
         0
    0.2508
    0.2478];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = ldesAvgTableRR.RespT(:);
Rsaved = [
         0
    0.5736
    0.5683];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ldesAvgTableRR.Tput(:);
Tsaved = [
    1.0000
    0.5031
    0.5031];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

lineTestAssertChecks(testName, gsFailures);
