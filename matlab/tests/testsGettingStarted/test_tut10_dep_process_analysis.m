evalc('tut10_dep_process_analysis');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end
SCVdEstsaved = 0.819332306743409;
SCVdsaved = 0.875000000000885;
TOL=1e-2; % 1-percent tolerance

try
    assert(max(max(abs(SCVdEst-SCVdEstsaved)))<TOL,sprintf('%s changed on SCVdEst.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(SCVdEst),mat2str(SCVdEstsaved))
end

try
    assert(max(max(abs(SCVd-SCVdsaved)))<TOL,sprintf('%s changed on SCVd.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(SCVd),mat2str(SCVdsaved))
end

lineTestAssertChecks(testName, gsFailures);
