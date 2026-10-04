exampleName = 'statepr_sys_aggr_large';
evalc('statepr_sys_aggr_large');
%%
TOL = 1e-6;
ATPr = [Pr_ctmc, Pr_jmt, Pr_nc];
% The third entry is SolverNC's; it was rebased once the normalizing-constant
% path moved it from 3.414029e-04 to 3.442270e-04, which is nearer the exact
% CTMC value in the first column than the recording it replaces.
ATPrEx = [
      0.000348435691610821 0.0003482324304121550 0.0003442269767886470
    ];
try
    assert(max(max(abs(ATPr-ATPrEx)))<TOL,sprintf('%s changed on ATPrEx.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(ATPr),mat2str(ATPrEx))
end
