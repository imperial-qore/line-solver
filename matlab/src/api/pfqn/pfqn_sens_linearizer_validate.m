function pfqn_sens_linearizer_validate()
%{
%{
 % @file pfqn_sens_linearizer_validate.m
 % @brief Validation harness for pfqn_sens_linearizer, the LINEARIZER-2 /
 %        LINEARIZER-3 moment approximation of Strelen (1990), Section 5.
 %
 %        This routine is an APPROXIMATION, so it must not be held to machine
 %        precision. The checks are therefore of three kinds:
 %
 %          A. accuracy bands against the exact pfqn_sens_mom, on models small
 %             enough for the exact lattice. The reference reports relative
 %             errors below 2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on E[Q^3]
 %             over its own 51 networks; the bands asserted here are of that
 %             order. This is the only meaningful statement of correctness for
 %             an approximation: it must track the exact answer, not equal it.
 %          B. exactness where the approximation degenerates. At a population of
 %             one job the CORE estimate of the queue lengths one job down is
 %             identically zero whatever the delta terms are, so the Linearizer
 %             equations coincide with the exact MVA and every moment must match
 %             pfqn_sens_mom to roundoff. This pins the derivative algebra
 %             independently of the heuristic.
 %          C. structural invariants that hold for any population: the mean
 %             queue lengths conserve the population, and the mean queue lengths
 %             agree with LINE's own pfqn_linearizer, which runs the same
 %             heuristic without derivatives.
 %
 %        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
 %        and its Linearizer", Performance Evaluation 11:127-142, 1990.
%}
%}
rng(7);
tolExact = 1e-9;    % B: the one-job case must be exact
tolPop   = 1e-8;    % C: population conservation
tolLin   = 5e-2;    % C: agreement with LINE's pfqn_linearizer on the means
% A: the bands are the accuracy the reference itself claims in Section 5 over
% its 51 networks, so this asserts that our port reproduces the paper's own
% accuracy statement rather than some slacker figure. rng is fixed, so the model
% set is deterministic and these are hard regression guards.
bandM    = 0.021;   % r(Q)   < 2.1% in the reference
bandM2   = 0.041;   % r(Q^2) < 4.1% in the reference
bandM3   = 0.062;   % r(Q^3) < 6.2% in the reference
% The reference does not report an error on the variance. It is naturally larger
% than the one on E[Q^2] because Var = E[Q^2] - E[Q]^2 is a difference of larger
% numbers, so relative error is amplified; banded here only to catch regressions.
bandVar  = 0.08;

errExact = 0; errPop = 0; errLin = 0;
worstM = 0; worstM2 = 0; worstM3 = 0; worstVar = 0;
nA = 0;

% =====================================================================
% A / C. random closed models against the exact moment analysis
% =====================================================================
for trial = 1:60
    M = randi([2 4]);
    R = randi([1 2]);
    L = 0.2 + rand(M,R);
    N = randi([1 4],1,R);
    if mod(trial,2) == 0
        Z = 0.5 + rand(1,R);
    else
        Z = zeros(1,R);
    end

    app = pfqn_sens_linearizer(L,N,Z);
    ex  = pfqn_sens_mom(L,N,Z);

    worstM  = max(worstM,  relerr(app.m,  ex.m));
    worstM2 = max(worstM2, relerr(app.M2, ex.M2));
    worstM3 = max(worstM3, relerr(app.M3, ex.M3));
    worstVar = max(worstVar, relerr(app.Var, ex.Var));
    nA = nA + 1;

    % ---- C. population conservation --------------------------------
    for r = 1:R
        inNet = sum(app.Q(:,r));
        inDelay = app.X(r) * Z(r);
        errPop = max(errPop, abs(inNet + inDelay - N(r)) / max(1,N(r)));
    end

    % ---- C. means against LINE's own Linearizer --------------------
    QL = pfqn_linearizer(L,N,Z,ones(M,1),1e-10,500);
    errLin = max(errLin, relerr(sum(app.Q,2), sum(QL,2)));
end

% =====================================================================
% B. one job: the Linearizer equations degenerate to the exact MVA
% =====================================================================
for trial = 1:12
    M = randi([2 4]);
    L = 0.2 + rand(M,1);
    Z = (mod(trial,2)==0) * (0.4 + rand);
    app = pfqn_sens_linearizer(L,1,Z);
    ex  = pfqn_sens_mom(L,1,Z);
    errExact = max([errExact, relerr(app.m,ex.m), relerr(app.Var,ex.Var), ...
                    relerr(app.M2,ex.M2), relerr(app.M3,ex.M3), ...
                    relerr(app.Cov,ex.Cov)]);
end

fprintf('\n=== pfqn_sens_linearizer validation ===\n');
fprintf('  A. accuracy vs exact pfqn_sens_mom (%d models), max relative error:\n', nA);
fprintf('       E[Q]   : %6.3f%%   (band %.1f%%)\n', 100*worstM,  100*bandM);
fprintf('       Var[Q] : %6.3f%%   (band %.1f%%, not reported by the paper)\n', 100*worstVar, 100*bandVar);
fprintf('       E[Q^2] : %6.3f%%   (band %.1f%%)\n', 100*worstM2, 100*bandM2);
fprintf('       E[Q^3] : %6.3f%%   (band %.1f%%)\n', 100*worstM3, 100*bandM3);
fprintf('  B. one job, exact vs pfqn_sens_mom      : %.3e  (tol %.1e)\n', errExact, tolExact);
fprintf('  C. population conservation              : %.3e  (tol %.1e)\n', errPop, tolPop);
fprintf('  C. means vs pfqn_linearizer             : %.3e  (tol %.1e)\n', errLin, tolLin);

ok = worstM <= bandM && worstM2 <= bandM2 && worstM3 <= bandM3 && ...
     worstVar <= bandVar && errExact <= tolExact && errPop <= tolPop && ...
     errLin <= tolLin;
if ~ok
    error('pfqn_sens_linearizer_validate:mismatch','one or more checks exceeded tolerance');
end
fprintf('  ALL CHECKS PASSED\n');
end

% =========================================================================
function e = relerr(a,b)
a = a(:); b = b(:);
d = abs(a-b);
scale = max(1, max(abs(a),abs(b)));
e = max(d ./ scale);
if isempty(e)
    e = 0;
end
end
