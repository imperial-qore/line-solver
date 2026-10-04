exampleName = 'ld_joint_dependence';
evalc('ld_joint_dependence');
TOL=1e-2; % 1-percent tolerance

%%
QN = jdAvgTableJD.QLen(:);
Qsaved = [0.889832793061531;0.528026484814507;15.1101672069385;7.47197351518549];
try
      assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = jdAvgTableJD.Util(:);
% Util rows 3-4 are the class-dependent multiserver station, normalized by the
% declared peak rate (Util = T*S/peak, peak=c=2), matching the T*S/c convention.
Usaved = [0.889832793063437;0.528026484815638;0.667374320414279;0.330016418819882];
try
      assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = jdAvgTableJD.RespT(:);
Rsaved = [0.999999999997858;1.99999999999572;16.9809062160078;28.3015103600129];
try
      assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

WN = jdAvgTableJD.ResidT(:);
Wsaved = [0.999999999997858;1.99999999999572;16.9809062160078;28.3015103600129];
try
      assert(max(max(abs(WN-Wsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(Wsaved))
end

TN = jdAvgTableJD.Tput(:);
Tsaved = [0.889832793063437;0.264013242407819;0.889832793063437;0.264013242407819];
try
      assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on Tput.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
%%
QN = jdAvgTableCTMC.QLen(:);
Qsaved = [
    0.892310275732868
    0.529227669119806
    15.107689724267134
    7.470772330880195
    ];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = jdAvgTableCTMC.Util(:);
Usaved = [
    0.892310275732868
    0.529227669119806
    0.669232706800043
    0.330767293199957
    ];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = jdAvgTableCTMC.RespT(:);
Rsaved =    [
    1.000000000000000
    2.000000000000000
    16.930982568647554
    28.232735235871434
    ];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

WN = jdAvgTableCTMC.ResidT(:);
Wsaved =    [
    1.000000000000000
    2.000000000000000
    16.930982568647554
    28.232735235871434
    ];
try
    assert(max(max(abs(WN-Wsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(Wsaved))
end

TN = jdAvgTableCTMC.Tput(:);
Tsaved = [
    0.892310275732868
    0.264613834559903
    0.892310275733391
    0.264613834559966
    ];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end