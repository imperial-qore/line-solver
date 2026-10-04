exampleName = 'ld_class_dependence';
evalc('ld_class_dependence');
TOL=1e-2; % 1-percent tolerance

%% product-form class dependence: beta_{i,r} reads only its own marginal,
%% class 1 scales up to c servers in n_{i,1}, class 2 is single-server.
QN = cdAvgTableCD.QLen(:);
Qsaved = [0.880423059243366;0.270071343739812;15.1195774807707;7.72992893234569];
try
      assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = cdAvgTableCD.Util(:);
% Util rows 3-4 are the class-dependent station, normalized by the declared
% per-class peak rate (Util = T*S/peak, peak=[c 1]).
Usaved = [0.880423029528273;0.270071334419442;0.660317272146205;0.337589168024303];
try
      assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = cdAvgTableCD.RespT(:);
% Rebased after the AMVA arrival queue dropped the arriving chain: MATLAB, the
% JAR and the C++ port all return this vector.
Rsaved = [1.00000003375724;2.00000006897116;17.1784625333749;57.2084638172119];
try
      assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

WN = cdAvgTableCD.ResidT(:);
Wsaved = [1.00000003375724;2.00000006897116;17.1784625333749;57.2084638172119];
try
      assert(max(max(abs(WN-Wsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(Wsaved))
end

TN = cdAvgTableCD.Tput(:);
Tsaved = [0.880423029528273;0.135035667209721;0.880423029528273;0.135035667209721];
try
      assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on Tput.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
%%
QN = cdAvgTableCTMC.QLen(:);
Qsaved = [0.881963178149607;0.270822093109951;15.1180368218504;7.72917790689005];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = cdAvgTableCTMC.Util(:);
Usaved = [0.881963178149607;0.270822093109951;0.661472383612429;0.338527616387571];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = cdAvgTableCTMC.RespT(:);
Rsaved = [1.000000000000000;2.000000000000000;17.14134693646;57.0793750105835];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

WN = cdAvgTableCTMC.ResidT(:);
Wsaved = [1.000000000000000;2.000000000000000;17.14134693646;57.0793750105835];
try
    assert(max(max(abs(WN-Wsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(Wsaved))
end

TN = cdAvgTableCTMC.Tput(:);
Tsaved = [0.881963178149607;0.135411046554976;0.881963178149905;0.135411046555028];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
