function x0 = qrf_noblo_start(nonlcon, n)
% QRF_NOBLO_START  A point of the QRF no-blocking polytope, for fmincon's start.
%
%   X0 = QRF_NOBLO_START(NONLCON, N) recovers the affine form of the
%   constraint callback NONLCON (which returns [C, CEQ] with C <= 0 and
%   CEQ == 0 for a decision vector of length N) and returns a point of the
%   polytope obtained by a phase-1 LP.
%
%   The all-zero vector that these methods used as a start violates ONE
%   (normalization) by a full unit and COR1 by N^2. fmincon then spends its
%   iteration budget restoring feasibility and terminates at a point that is
%   still infeasible: on a 2-station cyclic MAP/MAP network at N=2 it returned
%   a utilization ABOVE the LP maximum of the same polytope, and the queue
%   lengths did not conserve the population. Every constraint of this model is
%   linear, so a phase-1 LP lands on the polytope exactly and costs a fraction
%   of one fmincon iteration.
%
%   The phase 1 minimises 0.5*||x||^2 over the polytope rather than a null
%   objective, so the returned point is the unique minimum-norm feasible point
%   instead of an arbitrary vertex. (It is solved with QUADPROG, not LINPROG:
%   LINPROG errors with "Unrecognized field name optimstatus" on R2025a even
%   for a two-variable problem.)
%
%   The polytope is nonempty for every well-posed instance, so infeasibility
%   here is a modelling error and is raised rather than replaced by an
%   arbitrary point.

[c0, ceq0] = nonlcon(zeros(n,1));
c0 = full(c0(:));
ceq0 = full(ceq0(:));
A = zeros(numel(c0), n);
Aeq = zeros(numel(ceq0), n);
for col = 1:n
    e = zeros(n,1);
    e(col) = 1;
    [c1, ceq1] = nonlcon(e);
    A(:,col) = full(c1(:)) - c0;
    Aeq(:,col) = full(ceq1(:)) - ceq0;
end
b = -c0;
beq = -ceq0;

% Drop linearly dependent equality rows before the phase 1. The block is
% heavily redundant (SYMMETRY states every station pair twice, and the ZERO,
% MARGINALS and UEFF families overlap), carrying roughly twice as many rows as
% its rank, and quadprog stalls on the rank-deficient system: on a 2-station
% Erlang-2/Erlang-2 network at N=3 it stops at an equality residual of 2e-5
% instead of machine precision. Dropping dependent rows changes no feasible
% point, since they are exact linear combinations of the retained ones.
keep = qrf_independent_rows(Aeq);
AeqR = Aeq(keep,:);
beqR = beq(keep);
if rank([AeqR beqR]) > numel(keep)
    error('qrf_noblo_start:inconsistent', ...
          'QRF no-blocking equality system is inconsistent');
end

opts = optimset('Display','off','MaxIter',1000,'TolFun',1e-12,'TolCon',1e-12);
x0 = quadprog(speye(n), zeros(n,1), A, b, AeqR, beqR, ...
              zeros(n,1), ones(n,1), [], opts);

% Judge the phase 1 by the residual of the point it returned, not by the
% exit flag: quadprog reports flag 0 (iteration limit) on instances where it
% has nevertheless landed on the polytope to machine precision.
if isempty(x0)
    error('qrf_noblo_start:infeasible', ...
          'QRF no-blocking polytope phase 1 returned no point');
end
x0 = x0(:);
req = 0; rub = 0;
if ~isempty(Aeq), req = max(abs(Aeq*x0 - beq)); end
if ~isempty(A), rub = max(A*x0 - b); end
if req > 1e-8 || rub > 1e-8
    error('qrf_noblo_start:infeasible', ...
          ['QRF no-blocking polytope phase 1 did not reach feasibility ' ...
           '(max equality residual %.3e, max inequality residual %.3e)'], req, rub);
end
end
