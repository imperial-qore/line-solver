exampleName = 'oqn_vsinks';
evalc('oqn_vsinks');
TOL=1e-2; % 1-percent tolerance

QN = AvgTable.QLen(:);
Qsaved = [
                   0
                   0
   0.010204081106010
   0.010204081516070];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTable.Util(:);
Usaved = [
                   0
                   0
   0.010000000000000
   0.010000000000000];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgTable.RespT(:);
Rsaved = [
                   0
                   0
   0.010204081106010
   0.010204081516070];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = AvgTable.Tput(:);
Tsaved = [
     1
     1
     1
     1];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end

QN = AvgNodeTable.QLen(:);
Qsaved = [
                   0
                   0
   0.010204081106010
   0.010204081516070
                   0
                   0
                   0
                   0
                   0
                   0];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgNodeTable.Util(:);
Usaved = [
                   0
                   0
   0.010000000000000
   0.010000000000000
                   0
                   0
                   0
                   0
                   0
                   0];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgNodeTable.RespT(:);
Rsaved = [
                   0
                   0
   0.010204081106010
   0.010204081516070
                   0
                   0
                   0
                   0
                   0
                   0];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = AvgNodeTable.Tput(:);
Tsaved = [
   1.000000000000000
   1.000000000000000
   1.000000000000000
   1.000000000000000
                   0
                   0
   0.600000000000000
   0.100000000000000
   0.400000000000000
   0.900000000000000];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
%%
QN = AvgTableMAM.QLen(:);
Qsaved = [0;0;0.0102040816326531;0.0102040816326531];
try
      assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTableMAM.Util(:);
Usaved = [0;0;0.01;0.01];
try
      assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgTableMAM.RespT(:);
Rsaved = [0;0;0.0102040816326531;0.0102040816326531];
try
      assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

WN = AvgTableMAM.ResidT(:);
Wsaved = [0;0;0.0102040816326531;0.0102040816326531];
try
      assert(max(max(abs(WN-Wsaved)))<TOL,sprintf('%s changed on ResidT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(WN),mat2str(Wsaved))
end

TN = AvgTableMAM.Tput(:);
Tsaved = [1;1;1;1];
try
      assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on Tput.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
%%

QN = AvgNodeTableMAM.QLen(:);
Qsaved = [0;0;0.0102040816326531;0.0102040816326531;0;0;0;0;0;0];
try
      assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QLen.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgNodeTableMAM.Util(:);
Usaved = [0;0;0.01;0.01;0;0;0;0;0;0];
try
      assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on Util.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgNodeTableMAM.RespT(:);
Rsaved = [0;0;0.0102040816326531;0.0102040816326531;0;0;0;0;0;0];
try
      assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RespT.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = AvgNodeTableMAM.Tput(:);
Tsaved = [1;1;1;1;0;0;0.6;0.1;0.4;0.9];
try
      assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on Tput.',exampleName));
catch me
      fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end