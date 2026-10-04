% test_bolch_sum_ch10_12 - SUM/ESUM (closing method) vs Bolch Example 10.12.
% Open non-product-form network (Fig. 10.20, Tables 10.18/10.19): N=4
% single-server FCFS nodes, lambda_0=3, c2_0=1.5. Solved with the closing
% method + SUM (api/sum: sum_closing). Reference values are Bolch Table 10.19
% (T and Kbar vs the closed-population truncation K_closed). p41=0.1 is used
% (the book's printed 0.4 is a documented erratum that makes node 4 unstable;
% see line-stash bolch/ERRATUM.md #6).
exampleName = 'bolch_sum_ch10_12';

p41 = 0.1;
e  = [1; 3/7; 7/9; 1] / (1 - p41);   % visit ratios per external arrival
mu = [9; 10; 12; 4];
c2 = [0.5; 0.8; 2.4; 4.0];
L  = e ./ mu;

Kc_list = [100 200 500 1000 5000 10000];
T_book  = [3.94 4.16 4.29 4.34 4.37 4.37];   % Bolch Table 10.19 throughput
K_book  = [11.83 12.49 12.88 13.02 13.12 13.12]; % Bolch Table 10.19 mean population

% X must converge to the external arrival rate lambda_0=3; T/Kbar must fall
% within the documented ESUM accuracy band (~10%) of the book values.
XTOL  = 1e-2;
RELTOL = 0.10;
for idx = 1:numel(Kc_list)
    % 5th output is the mean system response time T (book Table 10.19 "T");
    % by Little's law T = Kbar/X. X (1st output) must reach lambda_0.
    [X, Q, ~, ~, T] = sum_closing(3, 1.5, L, [1;1;1;1], c2, Inf, 0, Kc_list(idx));
    Kbar = sum(Q);
    try
        assert(abs(X - 3) < XTOL, ...
            sprintf('%s: X did not reach lambda_0=3 at K=%d (got %.4f).', exampleName, Kc_list(idx), X));
        assert(abs(T - T_book(idx)) / T_book(idx) < RELTOL, ...
            sprintf('%s: T out of ESUM band at K=%d (LINE %.4f vs book %.4f).', exampleName, Kc_list(idx), T, T_book(idx)));
        assert(abs(Kbar - K_book(idx)) / K_book(idx) < RELTOL, ...
            sprintf('%s: Kbar out of ESUM band at K=%d (LINE %.4f vs book %.4f).', exampleName, Kc_list(idx), Kbar, K_book(idx)));
    catch me
        fprintf('Assertion failed in %s.\n%s\n', exampleName, me.message);
    end
end
