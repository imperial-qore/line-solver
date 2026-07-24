% Run all inference tests
% Requires lineStart to have been called first
tests = {'test_infer_ekf','test_infer_mcmc','test_infer_erps','test_infer_mle','test_infer_ubo','test_infer_ubr','test_infer_qmle','test_infer_mlps','test_infer_fmlps','test_infer_gibbs'};
for i=1:length(tests)
    fprintf('\n========== %s ==========\n', tests{i});
    tic;
    try
        out = evalc(tests{i});
        elapsed = toc;
        fprintf('PASS (%.1fs)\n', elapsed);
        % print all output lines containing estimates
        lines = strsplit(out, '\n');
        printNext = false;
        for l=1:length(lines)
            ln = strtrim(lines{l});
            if contains(ln, 'estVal') || contains(ln, 'Estimated')
                fprintf('%s\n', ln);
                printNext = true;
            elseif printNext && ~isempty(ln)
                fprintf('%s\n', ln);
                printNext = false;
            else
                printNext = false;
            end
        end
    catch e
        elapsed = toc;
        fprintf('FAIL (%.1fs): %s\n', elapsed, e.message);
    end
end
fprintf('\nDONE\n');
