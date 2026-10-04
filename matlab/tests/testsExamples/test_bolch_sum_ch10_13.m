% test_bolch_sum_ch10_13 - SUM (closing method) vs Bolch Example 10.13.
% Mixed non-product-form network (Figs. 10.21/10.22, Tables 10.20/10.21):
% N=5 multiserver nodes, R=2 classes (class 1 closed K_1=9; class 2 open
% lambda_01=5, c2_01=0.7). Solved with the closing method + SUM (api/sum:
% sum_closing). Reference values are Bolch Table 10.21 (closing method CLS):
% rho and Kbar per node, reproduced at printed precision.
exampleName = 'bolch_sum_ch10_13';

% visit ratios (class 1 normalized at node 1; class 2 per external arrival)
e1 = [1; 10/7; 10/7; 3/7; 3/7];
e2 = [1; 1.25; 10/7; 0.25; 3/7];
mu = [4; 4; 6; 4; 5];
m  = [3; 4; 3; 2; 2];
c2 = [0.4; 0.3; 0.3; 0.4; 0.5];
L  = [e1 ./ mu, e2 ./ mu];

rho_book = [0.84 0.85 0.80 0.43 0.43];   % Bolch Table 10.21 utilization
K_book   = [5.2 5.8 4.1 1.0 1.0];        % Bolch Table 10.21 mean population

[X, Q, U] = sum_closing([0 5], [1 0.7], L, m, [c2 c2], [9 Inf], [0 0], 500);
rho  = sum(U, 2)';
Kbar = sum(Q, 2)';

RHOTOL = 5e-3;   % printed to 2 decimals in the book
KTOL   = 5e-2;
try
    assert(abs(X(2) - 5) < 1e-2, ...
        sprintf('%s: open class throughput did not reach lambda_01=5 (got %.4f).', exampleName, X(2)));
    for i = 1:5
        assert(abs(rho(i) - rho_book(i)) < RHOTOL, ...
            sprintf('%s: rho(node %d) LINE %.4f vs book %.2f.', exampleName, i, rho(i), rho_book(i)));
        assert(abs(Kbar(i) - K_book(i)) < KTOL, ...
            sprintf('%s: Kbar(node %d) LINE %.4f vs book %.1f.', exampleName, i, Kbar(i), K_book(i)));
    end
catch me
    fprintf('Assertion failed in %s.\n%s\n', exampleName, me.message);
end
