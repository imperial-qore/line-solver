exampleName = 'statepr_aggr_large';
evalc('statepr_aggr_large');
TOL = 1e-6;
Pr_ctmcsaved = 0.005511140445099;
try
    assert(max(max(abs(Pr_ctmc-Pr_ctmcsaved)))<TOL,sprintf('%s changed on Pr_ctmc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pr_ctmc),mat2str(Pr_ctmcsaved))
end
Pr_ncsaved = 0.005511140445099;
try
    assert(max(max(abs(Pr_nc-Pr_ncsaved)))<TOL,sprintf('%s changed on Pr_nc.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(Pr_nc),mat2str(Pr_ncsaved))
end
