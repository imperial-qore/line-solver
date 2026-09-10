%{ @file dph_from_dist.m
 %  @brief Exact discrete phase-type representation of a lattice-valued law
 %
 %  @author LINE Development Team
%}

%{
 % @brief Exact discrete phase-type representation of a lattice-valued law
 %
 % @details
 % Returns (alpha, A) with P[X=k] = alpha*A^(k-1)*a, a = e - A*e, k = 1,2,...,
 % measuring X in slots. The three families below are represented EXACTLY, not
 % moment-matched: Geometric by its single-phase chain, Det(k) by the k-step
 % chain, DiscreteUniform by the hazard chain h_j = 1/(hi-j+1) on j >= lo. A
 % family outside this set is an error rather than a fit, because a fitted
 % surrogate would silently leave the lattice the caller is relying on.
 %
 % @par Syntax:
 % @code
 % [alpha, A] = dph_from_dist(procType, meanSlots, scv)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>procType<td>ProcessType id (GEOMETRIC, DET, DUNIFORM)
 % <tr><td>meanSlots<td>Mean of the law expressed in slots
 % <tr><td>scv<td>Squared coefficient of variation (used by DUNIFORM only)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>alpha<td>1xm initial phase probability row vector
 % <tr><td>A<td>mxm substochastic transient matrix
 % </table>
%}
function [alpha, A] = dph_from_dist(procType, meanSlots, scv)

tol = GlobalConstants.FineTol;

switch procType
    case ProcessType.GEOMETRIC
        p = 1 / meanSlots;
        if p > 1 + tol || p <= 0
            line_error(mfilename, 'Geometric with mean %g slots is outside the support {1,2,...}.', meanSlots);
        end
        p = min(1, p);
        alpha = 1;
        A = 1 - p;
    case ProcessType.DET
        k = round(meanSlots);
        if abs(meanSlots - k) > tol * max(1, meanSlots) || k < 1
            line_error(mfilename, 'Det of %g slots is not a positive integral number of slots.', meanSlots);
        end
        alpha = zeros(1, k); alpha(1) = 1;
        A = zeros(k, k);
        for i = 1:(k-1)
            A(i, i+1) = 1;
        end
    case ProcessType.DUNIFORM
        varSlots = scv * meanSlots^2;
        width = sqrt(max(0, 12 * varSlots + 1)) - 1;
        lo = round(meanSlots - width / 2);
        hi = round(meanSlots + width / 2);
        if lo < 1 || hi < lo
            line_error(mfilename, 'DiscreteUniform spanning [%d,%d] slots is outside the support {1,2,...}.', lo, hi);
        end
        alpha = zeros(1, hi); alpha(1) = 1;
        A = zeros(hi, hi);
        for j = 1:(hi-1)
            if j < lo
                h = 0;
            else
                h = 1 / (hi - j + 1);
            end
            A(j, j+1) = 1 - h;
        end
        % the last phase always absorbs, so row hi stays zero
    otherwise
        line_error(mfilename, ['ProcessType %s has no exact discrete phase-type representation. ' ...
            'The discrete-time path accepts Geometric, Det on the slot lattice, ' ...
            'DiscreteUniform and DMAP.'], ProcessType.toText(procType));
end

end
