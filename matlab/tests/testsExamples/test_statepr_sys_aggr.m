exampleName = 'statepr_sys_aggr';
evalc('statepr_sys_aggr');
TOL = 1e-6;
Pr_ctmcsaved = 2.626044433948373e-04;
try
    assert(max(max(abs(Pr_ctmc-Pr_ctmcsaved)))<TOL,sprintf('%s changed on Pr_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pr_ctmc),mat2str(Pr_ctmcsaved))
end
Pr_ncsaved = 2.626044433948373e-04;
try
    assert(max(max(abs(Pr_nc-Pr_ncsaved)))<TOL,sprintf('%s changed on Pr_nc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pr_nc),mat2str(Pr_ncsaved))
end
