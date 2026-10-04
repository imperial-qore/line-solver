function residual = qrf_lp_residual(x, Aeq, beq, Aineq, bineq, lb, ub)
% QRF_LP_RESIDUAL Maximum constraint violation of a point of a QRF bound LP.
%
%   residual = QRF_LP_RESIDUAL(x, Aeq, beq, Aineq, bineq, lb, ub) returns the
%   largest violation of the equality rows, the inequality rows and the box, as
%   an absolute quantity. A point returned by linprog with a nonzero exitflag
%   can be finite yet infeasible; qrf_bas and qrf_rsrd use this to tell a
%   solution of the LP from a point that merely came back from the solver.
%
%   See also QRF_BAS, QRF_RSRD.

residual = 0;
if ~isempty(Aeq)
    residual = max(residual, max(abs(Aeq * x - beq(:))));
end
if ~isempty(Aineq)
    residual = max(residual, max(Aineq * x - bineq(:)));
end
if ~isempty(lb)
    residual = max(residual, max(lb(:) - x));
end
if ~isempty(ub)
    finiteUb = isfinite(ub(:));
    if any(finiteUb)
        residual = max(residual, max(x(finiteUb) - ub(finiteUb)));
    end
end
residual = max(residual, 0);
end
