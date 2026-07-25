%{ @file sn_is_phasetype.m
 %  @brief Tests whether a process representation admits a phase-type reading
 %
 %  @author LINE Development Team
%}

%{
 % @brief Tests whether {D0,D1,...} is a Markovian (phase-type / MAP) pair
 %
 % @details
 % A representation is Markovian when D0 has nonnegative off-diagonal entries,
 % every D_k with k >= 1 is nonnegative, and the entry vector pie is
 % nonnegative. Exactly under those conditions do sn.mu, sn.phi and sn.pie carry
 % their probabilistic reading (mu_i = -D0(i,i) is a rate, phi_i is a completion
 % probability, pie is a distribution over phases), which is what the CTMC state
 % space, SSA and the fluid ODEs consume.
 %
 % A matrix-exponential (ME) or rational (RAP) process fails the test: its
 % moments, transforms and aggregated stationary measures remain exact, but the
 % per-phase quantities are signed. See _kb/04-networkstruct.md.
 %
 % @par Syntax:
 % @code
 % tf = sn_is_phasetype(MAP)
 % tf = sn_is_phasetype(MAP, pie)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAP<td>Cell array {D0,D1,...} of square matrices of equal size
 % <tr><td>pie<td>Optional entry vector to test for nonnegativity
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>tf<td>true when the representation is Markovian, false otherwise
 % </table>
%}
function tf = sn_is_phasetype(MAP, pie)

tol = GlobalConstants.Zero;

tf = true;

% An empty, scalar-parameter or NaN-carrying entry describes a disabled or
% not-yet-Markovian process; there is no phase decomposition to invalidate, so
% the caller is not blocked by this test.
if isempty(MAP) || ~iscell(MAP) || numel(MAP) < 2
    return;
end

D0 = MAP{1};
if ~isnumeric(D0) || isempty(D0) || any(isnan(D0(:)))
    return;
end

n = size(D0,1);
if size(D0,2) ~= n
    return;
end

% Off-diagonal entries of D0 are transition rates between phases.
offdiag = D0 - diag(diag(D0));
if any(offdiag(:) < -tol)
    tf = false;
    return;
end

% D1 and any further D_k are jump matrices and must be nonnegative.
for k = 2:numel(MAP)
    Dk = MAP{k};
    if ~isnumeric(Dk) || any(isnan(Dk(:)))
        continue;
    end
    if any(Dk(:) < -tol)
        tf = false;
        return;
    end
end

if nargin > 1 && ~isempty(pie) && isnumeric(pie) && ~any(isnan(pie(:)))
    if any(pie(:) < -tol)
        tf = false;
        return;
    end
end
end
