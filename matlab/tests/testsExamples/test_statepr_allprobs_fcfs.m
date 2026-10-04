exampleName = 'statepr_allprobs_fcfs';
evalc('statepr_allprobs_fcfs');
%%
TOL = 1e-6;
Pmarga_ctmcsaved = 0.307692307692308;
try
    assert(max(max(abs(Pmarga_ctmc-Pmarga_ctmcsaved)))<TOL,sprintf('%s changed on Pmarga_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarga_ctmc),mat2str(Pmarga_ctmcsaved))
end
Pmarga_ncsaved = 0.307692307692308;
try
    assert(max(max(abs(Pmarga_nc-Pmarga_ncsaved)))<TOL,sprintf('%s changed on Pmarga_nc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarga_nc),mat2str(Pmarga_ncsaved))
end
Pmarga_ssasaved = 0.308420743700602;
try
    assert(max(max(abs(Pmarga_ssa-Pmarga_ssasaved)))<TOL,sprintf('%s changed on Pmarga_ssa.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarga_ssa),mat2str(Pmarga_ssasaved))
end
Pmarga_jmtsaved = 0.308371392687713;
try
    assert(max(max(abs(Pmarga_jmt-Pmarga_jmtsaved)))<TOL,sprintf('%s changed on Pmarga_jmt.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarga_jmt),mat2str(Pmarga_jmtsaved))
end
Pmarg_ctmcsaved = 0.307692307692308;
try
    assert(max(max(abs(Pmarg_ctmc-Pmarg_ctmcsaved)))<TOL,sprintf('%s changed on Pmarg_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarg_ctmc),mat2str(Pmarg_ctmcsaved))
end
Pmarg_ssasaved = 0.308420743700602;
try
    assert(max(max(abs(Pmarg_ssa-Pmarg_ssasaved)))<TOL,sprintf('%s changed on Pmarg_ssa.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarg_ssa),mat2str(Pmarg_ssasaved))
end
Pjointa_ctmcsaved = 0.076923076923077;
try
    assert(max(max(abs(Pjointa_ctmc-Pjointa_ctmcsaved)))<TOL,sprintf('%s changed on Pjointa_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjointa_ctmc),mat2str(Pjointa_ctmcsaved))
end
Pjointa_ncsaved = 0.076190930869983;
try
    assert(max(max(abs(Pjointa_nc-Pjointa_ncsaved)))<TOL,sprintf('%s changed on Pjointa_nc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjointa_nc),mat2str(Pjointa_ncsaved))
end
Pjointa_ssasaved = 0.077739789086896;
try
    assert(max(max(abs(Pjointa_ssa-Pjointa_ssasaved)))<TOL,sprintf('%s changed on Pjointa_ssa.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjointa_ssa),mat2str(Pjointa_ssasaved))
end
Pjointa_jmtsaved = 0.076127055280772;
try
    assert(max(max(abs(Pjointa_jmt-Pjointa_jmtsaved)))<TOL,sprintf('%s changed on Pjointa_jmt.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjointa_jmt),mat2str(Pjointa_jmtsaved))
end
Pjoint_ctmcsaved = 0.076923076923077;
try
    assert(max(max(abs(Pjoint_ctmc-Pjoint_ctmcsaved)))<TOL,sprintf('%s changed on Pjoint_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjoint_ctmc),mat2str(Pjoint_ctmcsaved))
end
Pjoint_ssasaved = 0.077739789086896;
try
    assert(max(max(abs(Pjoint_ssa-Pjoint_ssasaved)))<TOL,sprintf('%s changed on Pjoint_ssa.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjoint_ssa),mat2str(Pjoint_ssasaved))
end