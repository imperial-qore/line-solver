exampleName = 'statepr_aggr';
evalc('statepr_aggr');
TOL = 1e-6;
Pmarg_ctmcsaved = 0.340000000000000;
try
    assert(max(max(abs(Pmarg_ctmc-Pmarg_ctmcsaved)))<TOL,sprintf('%s changed on Pmarg_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarg_ctmc),mat2str(Pmarg_ctmcsaved))
end
Pmarg_ncsaved = 0.340000000000000;
try
    assert(max(max(abs(Pmarg_nc-Pmarg_ncsaved)))<TOL,sprintf('%s changed on Pmarg_nc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pmarg_nc),mat2str(Pmarg_ncsaved))
end
