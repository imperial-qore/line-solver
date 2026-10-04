exampleName = 'oqn_trace_driven';
evalc('oqn_trace_driven');
TOL=1e-2; % 1-percent tolerance

%% Test AvgTable{1} (JMT)
QN = AvgTable{1}.QLen(:);
Qsaved = [
                  0
   0.112795903042015];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN (JMT).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTable{1}.Util(:);
Usaved = [
                   0
   0.102842847883759];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN (JMT).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgTable{1}.RespT(:);
Rsaved = [
                   0
   0.112746401331119];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN (JMT).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = AvgTable{1}.Tput(:);
Tsaved = [
   1.020491254100682
   1.021117709190517];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN (JMT).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end

%% Test AvgTable{2} (DES)
QN = AvgTable{2}.QLen(:);
Qsaved = [
                  0
   0.11483];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN (DES).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = AvgTable{2}.Util(:);
Usaved = [
                   0
   0.10127];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN (DES).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = AvgTable{2}.RespT(:);
Rsaved = [
                   0
   0.11215];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN (DES).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = AvgTable{2}.Tput(:);
Tsaved = [
   1
   1.002];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN (DES).',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
