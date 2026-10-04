evalc('tut02_mg1_multiclass_solvers');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end
TOL = 1e-2; % 1-percent tolerance

% ldesAvgTable
QN = ldesAvgTable.QLen(:);
Qsaved = [
    0
    0
    0.910773679443183
    0.429793424131847];
try
    assert(max(max(abs(QN - Qsaved))) < TOL, sprintf('%s changed on QN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(QN), mat2str(Qsaved))
end

UN = ldesAvgTable.Util(:);
Usaved = [
    0
    0
    0.500399676445432
    0.0503253211805998];
try
    assert(max(max(abs(UN - Usaved))) < TOL, sprintf('%s changed on UN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(UN), mat2str(Usaved))
end

RN = ldesAvgTable.RespT(:);
Rsaved = [
    0
    0
    1.83565283234187
    0.878411719852989];
try
    assert(max(max(abs(RN - Rsaved))) < TOL, sprintf('%s changed on RN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(RN), mat2str(Rsaved))
end

TN = ldesAvgTable.Tput(:);
Tsaved = [
    0.5
    0.5
    0.50053907319656
    0.493935919322689];
try
    assert(max(max(abs(TN - Tsaved))) < TOL, sprintf('%s changed on TN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(TN), mat2str(Tsaved))
end

% ctmcAvgTable2
QN = ctmcAvgTable2.QLen(:);
Qsaved = [
    0
    0
    0.567339369711028
    0.244555866368181];
try
    assert(max(max(abs(QN - Qsaved))) < TOL, sprintf('%s changed on QN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(QN), mat2str(Qsaved))
end

UN = ctmcAvgTable2.Util(:);
Usaved = [
    0
    0
    0.441071671607959
    0.0480923780382475];
try
    assert(max(max(abs(UN - Usaved))) < TOL, sprintf('%s changed on UN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(UN), mat2str(Usaved))
end

RN = ctmcAvgTable2.RespT(:);
Rsaved = [
    0
    0
    1.28627478532627
    0.513959982291845];
try
    assert(max(max(abs(RN - Rsaved))) < TOL, sprintf('%s changed on RN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(RN), mat2str(Rsaved))
end

TN = ctmcAvgTable2.Tput(:);
Tsaved = [
    0.441071671607959
    0.475826668990182
    0.441071671607959
    0.475826668990182];
try
    assert(max(max(abs(TN - Tsaved))) < TOL, sprintf('%s changed on TN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(TN), mat2str(Tsaved))
end

% ctmcAvgTable4
QN = ctmcAvgTable4.QLen(:);
Qsaved = [
    0
    0
    0.795796030011908
    0.375579325562571];
try
    assert(max(max(abs(QN - Qsaved))) < TOL, sprintf('%s changed on QN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(QN), mat2str(Qsaved))
end

UN = ctmcAvgTable4.Util(:);
Usaved = [
    0
    0
    0.491617843729365
    0.0501036818869567];
try
    assert(max(max(abs(UN - Usaved))) < TOL, sprintf('%s changed on UN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(UN), mat2str(Usaved))
end

RN = ctmcAvgTable4.RespT(:);
Rsaved = [
    0
    0
    1.61872893785766
    0.757634092533136];
try
    assert(max(max(abs(RN - Rsaved))) < TOL, sprintf('%s changed on RN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(RN), mat2str(Rsaved))
end

TN = ctmcAvgTable4.Tput(:);
Tsaved = [
    0.491617843729365
    0.495726537736479
    0.491617843729365
    0.495726537736479];
try
    assert(max(max(abs(TN - Tsaved))) < TOL, sprintf('%s changed on TN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(TN), mat2str(Tsaved))
end

% mamAvgTable
QN = mamAvgTable.QLen(:);
Qsaved = [
    0
    0
    0.876460023389536
    0.426995629286018];
try
    assert(max(max(abs(QN - Qsaved))) < TOL, sprintf('%s changed on QN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(QN), mat2str(Qsaved))
end

UN = mamAvgTable.Util(:);
Usaved = [
    0
    0
    0.5
    0.0505356058964822];
try
    assert(max(max(abs(UN - Usaved))) < TOL, sprintf('%s changed on UN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(UN), mat2str(Usaved))
end

RN = mamAvgTable.RespT(:);
Rsaved = [
    0
    0
    1.75292004677907
    0.853991258572036];
try
    assert(max(max(abs(RN - Rsaved))) < TOL, sprintf('%s changed on RN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(RN), mat2str(Rsaved))
end

TN = mamAvgTable.Tput(:);
Tsaved = [
    0.5
    0.5
    0.5
    0.5];
try
    assert(max(max(abs(TN - Tsaved))) < TOL, sprintf('%s changed on TN.', testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n', testName, me.message, mat2str(TN), mat2str(Tsaved))
end

lineTestAssertChecks(testName, gsFailures);
