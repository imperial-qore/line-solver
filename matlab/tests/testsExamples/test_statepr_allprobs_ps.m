    exampleName = 'statepr_allprobs_ps';
evalc('statepr_allprobs_ps');
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
Pmarga_jmtsaved = 0.307702307908087;
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
Pmarg_ssasaved = 0.308875397420987;
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
Pjointa_ncsaved = 0.07619047619047617;
try
    assert(max(max(abs(Pjointa_nc-Pjointa_ncsaved)))<TOL,sprintf('%s changed on Pjointa_nc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pjointa_nc),mat2str(Pjointa_ncsaved))
end
Pjointa_ssasaved = 0.0802855830031828;
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
