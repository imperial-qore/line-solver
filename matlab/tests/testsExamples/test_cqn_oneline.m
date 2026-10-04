exampleName = 'cqn_oneline';
evalc('cqn_oneline');
TOL=1e-2; % 1-percent tolerance

QN = avgTable.QLen(:);
Qsaved = [
    0.455045157326239
    0.915250653208990
    0.465046149794948
    0.935147406539621
    0.052560943314077
    0.053600791182090
    0.027347749564735
    0.096001149069299];
try
    assert(max(max(abs(QN-Qsaved)))<TOL,sprintf('%s changed on QN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(QN),mat2str(Qsaved))
end

UN = avgTable.Util(:);
Usaved = [
    0.455045157326239
    0.915250653208990
    0.465046149794948
    0.935147406539621
    0.050004962343543
    0.049741883326576
    0.025002481171771
    0.089535389987836];
try
    assert(max(max(abs(UN-Usaved)))<TOL,sprintf('%s changed on UN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(UN),mat2str(Usaved))
end

RN = avgTable.RespT(:);
Rsaved =    [
    91.000000000000000
    92.000000000000000
    93.000000000000000
    94.000000000000000
    10.511145464519025
    5.387893219701717
    5.469007131102578
    9.649931069056272];
try
    assert(max(max(abs(RN-Rsaved)))<TOL,sprintf('%s changed on RN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(RN),mat2str(Rsaved))
end

TN = avgTable.Tput(:);
Tsaved = [
    0.005000496234354
    0.009948376665315
    0.005000496234354
    0.009948376665315
    0.005000496234354
    0.009948376665315
    0.005000496234354
    0.009948376665315];
try
    assert(max(max(abs(TN-Tsaved)))<TOL,sprintf('%s changed on TN.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(TN),mat2str(Tsaved))
end
