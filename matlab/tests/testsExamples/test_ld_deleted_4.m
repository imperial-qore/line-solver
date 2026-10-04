exampleName = 'ld_heterogeneous_fcfs';
evalc('ld_heterogeneous_fcfs');
TOL=1e-2; % 1-percent tolerance

new_lldAvgTableNC = lldAvgTableNC;
new_lldAvgTableRD = lldAvgTableRD;
new_lldAvgTableNRL = lldAvgTableNRL;
new_lldAvgTableNRP = lldAvgTableNRP;
new_lldAvgTableMVALD = lldAvgTableMVALD;
new_lldAvgTableQD = lldAvgTableQD;
if exist('cdAvgTableCD','var')
    new_cdAvgTableCD = cdAvgTableCD;
end
new_lldAvgTableCTMC = lldAvgTableCTMC;
load(fullfile(fileparts(which('runTestExample')), 'testsExamples', [exampleName,'.mat']));


%%
QN = new_lldAvgTableNC.QLen(:);
QNsaved = lldAvgTableNC.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableNC.Util(:);
UNsaved = lldAvgTableNC.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableNC.RespT(:);
RNsaved = lldAvgTableNC.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableNC.ResidT(:);
WNsaved = lldAvgTableNC.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableNC.Tput(:);
TNsaved = lldAvgTableNC.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end

%%
QN = new_lldAvgTableRD.QLen(:);
QNsaved = lldAvgTableRD.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableRD.Util(:);
UNsaved = lldAvgTableRD.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableRD.RespT(:);
RNsaved = lldAvgTableRD.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableRD.ResidT(:);
WNsaved =  lldAvgTableRD.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableRD.Tput(:);
TNsaved = lldAvgTableRD.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end

%%
QN = new_lldAvgTableNRL.QLen(:);
QNsaved = lldAvgTableNRL.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableNRL.Util(:);
UNsaved = lldAvgTableNRL.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableNRL.RespT(:);
RNsaved = lldAvgTableNRL.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableNRL.ResidT(:);
WNsaved = lldAvgTableNRL.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableNRL.Tput(:);
TNsaved = lldAvgTableNRL.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end

%%
QN = new_lldAvgTableNRP.QLen(:);
QNsaved = lldAvgTableNRP.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableNRP.Util(:);
UNsaved = lldAvgTableNRP.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableNRP.RespT(:);
RNsaved = lldAvgTableNRP.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableNRP.ResidT(:);
WNsaved = lldAvgTableNRP.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableNRP.Tput(:);
TNsaved = lldAvgTableNRP.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end

%%
QN = new_lldAvgTableMVALD.QLen(:);
QNsaved = lldAvgTableMVALD.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableMVALD.Util(:);
UNsaved = lldAvgTableMVALD.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableMVALD.RespT(:);
RNsaved = lldAvgTableMVALD.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableMVALD.ResidT(:);
WNsaved = lldAvgTableMVALD.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableMVALD.Tput(:);
TNsaved = lldAvgTableMVALD.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end

%%
QN = new_lldAvgTableQD.QLen(:);
QNsaved = lldAvgTableQD.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableQD.Util(:);
UNsaved = lldAvgTableQD.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableQD.RespT(:);
RNsaved = lldAvgTableQD.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableQD.ResidT(:);
WNsaved = lldAvgTableQD.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableQD.Tput(:);
TNsaved = lldAvgTableQD.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on Tput.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end
%%
QN = new_lldAvgTableCTMC.QLen(:);
QNsaved = lldAvgTableCTMC.QLen(:);
try
    assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
end

UN = new_lldAvgTableCTMC.Util(:);
UNsaved = lldAvgTableCTMC.Util(:);
try
    assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
end

RN = new_lldAvgTableCTMC.RespT(:);
RNsaved = lldAvgTableCTMC.RespT(:);
try
    assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
end

WN = new_lldAvgTableCTMC.ResidT(:);
WNsaved = lldAvgTableCTMC.ResidT(:);
try
    assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
end

TN = new_lldAvgTableCTMC.Tput(:);
TNsaved = lldAvgTableCTMC.Tput(:);
try
    assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
end
%%
if exist('cdAvgTableCD','var')
    QN = new_cdAvgTableCD.QLen(:);
    QNsaved = cdAvgTableCD.QLen(:);
    try
        assert(max(max(abs(QN-QNsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
    catch me
        fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(QNsaved))
    end

    UN = new_cdAvgTableCD.Util(:);
    UNsaved = cdAvgTableCD.Util(:);
    try
        assert(max(max(abs(UN-UNsaved)))<TOL,sprintf('%s changed on Util.',exampleName));
    catch me
        fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(UNsaved))
    end

    RN = new_cdAvgTableCD.RespT(:);
    RNsaved = cdAvgTableCD.RespT(:);
    try
        assert(max(max(abs(RN-RNsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
    catch me
        fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(RNsaved))
    end

    WN = new_cdAvgTableCD.ResidT(:);
    WNsaved = cdAvgTableCD.ResidT(:);
    try
        assert(max(max(abs(WN-WNsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
    catch me
        fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(WNsaved))
    end

    TN = new_cdAvgTableCD.Tput(:);
    TNsaved = cdAvgTableCD.Tput(:);
    try
        assert(max(max(abs(TN-TNsaved)))<TOL,sprintf('%s changed on Tput.',exampleName));
    catch me
        fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(TNsaved))
    end
end