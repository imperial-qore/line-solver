[~,~,p_opt] = evalc('tut09_opt_load_balancing');
gsFailures = {}; % each catch below records its mismatch; lineTestAssertChecks fails the test at the end
p_opt_saved = 0.6105;
TOL=1e-2; % 1-percent tolerance
try
    assert(max(max(abs(p_opt-p_opt_saved)))<TOL,sprintf('%s changed on p_opt.',testName));
catch me
    gsFailures{end+1} = me.message;
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',testName,me.message,mat2str(p_opt),mat2str(p_opt_saved))
end

lineTestAssertChecks(testName, gsFailures);
