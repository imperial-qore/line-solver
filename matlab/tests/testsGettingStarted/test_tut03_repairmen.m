evalc('tut03_repairmen');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end

TOL=1e-2; % 1-percent tolerance

QN = ctmcAvgTable.QLen(:);
Qsaved = [
    2.6648
    0.3352];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = ctmcAvgTable.Util(:);
Usaved = [
    2.6648
    0.1666];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = ctmcAvgTable.RespT(:);
Rsaved = [
    2.0000
    0.2515];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = ctmcAvgTable.Tput(:);
Tsaved = [
    1.3324
    1.3324];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(TN),mat2str(Tsaved))
end

Qsaved = [   -8.0000    8.0000         0         0
    0.5000   -8.5000    8.0000         0
         0    1.0000   -5.0000    4.0000
         0         0    1.5000   -1.5000];
try
    assert(max(max(abs(InfGen-Qsaved)))<TOL,sprintf('%s changed on InfGen.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(InfGen),mat2str(Qsaved))
end

StateSpacesaved = [ 
     0     1     2
     1     0     2
     2     0     1
     3     0     0];
try
    assert(max(max(abs(StateSpace-StateSpacesaved)))<TOL,sprintf('%s changed on StateSpace.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(StateSpace),mat2str(StateSpacesaved))
end

lineTestAssertChecks(testName, gsFailures);
