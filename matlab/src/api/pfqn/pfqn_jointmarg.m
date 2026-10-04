%{
%{
 % @file pfqn_jointmarg.m
 % @brief Joint probability of the per-station TOTAL queue lengths.
%}
%}

%{
%{
 % @brief Joint probability of the per-station TOTAL queue lengths.
 % @fn pfqn_jointmarg(n, L, N, infset, lGN, engine)
 % @param n Per-station total queue lengths.
 % @param L Service demand matrix, infinite servers included as rows.
 % @param N Population vector.
 % @param infset Rows of L that are infinite-server stations (optional).
 % @param lGN Log normalizing constant at N (optional).
 % @param engine Permanent engine (optional).
 % @return pjoint Joint probability of the total queue lengths.
 % @return lpjoint Its logarithm.
%}
%}
function [pjoint, lpjoint] = pfqn_jointmarg(n, L, N, infset, lGN, engine)
% [PJOINT, LPJOINT] = PFQN_JOINTMARG(N, L, NPOP, INFSET, LGN, ENGINE)
%
% Joint probability that station I holds N(I) jobs IN TOTAL, all classes
% summed out, in a closed multiclass product-form network:
%
%   P(n_1,...,n_M) = perm(A) / ( prod_r N_r! * prod_{j in INFSET} n_j! * G(N) )
%
% with A the demand matrix whose column R is repeated N_r times and whose row
% I is repeated n_i times, so A is square of order sum(N). Unlike
% PFQN_JOINT, which takes the delay as a single aggregated row, every
% infinite-server station keeps its own row here and contributes its own
% 1/n_j! -- the queueing stations contribute the n_i! that the permanent
% identity supplies, the infinite servers do not.
%
% Input:
%   N       - (M x 1) per-station total queue lengths, infinite servers
%             included; sum(N) must equal sum(NPOP)
%   L       - (M x R) demand matrix, infinite-server rows included
%   NPOP    - (1 x R) per-class populations
%   INFSET  - indices of the rows of L that are infinite-server stations,
%             empty by default (every station is a queue)
%   LGN     - log normalizing constant; computed with PFQN_CA when omitted,
%             aggregating the infinite-server rows into the think time (which
%             is exact: the delay stations aggregate by the multinomial
%             theorem, so G does not depend on how they are split)
%   ENGINE  - 'exact' (default), 'spm', 'bethe', 'heur', 'huberlaw' or
%             'adapart'
%
% Output:
%   PJOINT  - joint probability
%   LPJOINT - its logarithm, which survives populations the probability
%             itself underflows at
%
% ENGINE 'spm' IS THE ONLY ONE THAT DOES NOT EXPAND THE MATRIX. The other
% approximate engines are handed A itself, of order SUM(N), where every column
% has multiplicity one and their accuracy decays with the population. PERM_SPM
% takes the ROW-replicated matrix with the class populations as its column
% multiplicities, which is the regime its expansion is asymptotically exact
% in: the Laplace integral has dimension R-1 whatever the population is, so
% the cost is independent of SUM(N) and the relative error is O((R-1)/min(N)).
% Measured on a 3-station 2-class model: 12.8% at N = (1,1), 4.2% at (3,3),
% 2.1% at (6,6), tracking 1/(8 min N_r). It degrades the other way round, when
% the CLASS COUNT grows at fixed population (2.7% at R = 2, 21% at R = 7, all
% at N_r = 3), because R-1 is the dimension being expanded in.
%
% The bias is nearly constant across the lattice, so a caller sweeping the
% whole state space and RENORMALIZING to sum to one keeps far less of it: the
% total variation distance against the exact law is 5.0e-3 at N = (1,1),
% 8.4e-4 at (3,3) and 4.3e-4 at (5,5), better than 'bethe' and 'heur' at every
% population measured. Cost against 'exact', which collapses the repeated rows
% and columns and so grows with M and R rather than with SUM(N): comparable at
% R = 2, 12x faster at R = 4 and 16x at R = 7.
%
% The identity holds for load-independent single-server queues plus infinite
% servers. Multiserver and load-dependent stations break the n_i! factor and
% are the caller's responsibility to exclude.
%
% ZERO ELEMENTS are safe under the exact engine and only under it: a station
% holding no jobs contributes no row, a class with no jobs contributes no
% column, a zero demand is an ordinary zero entry of A, and the permanent of
% the empty matrix is 1. The approximate engines are REFUSED on a matrix with
% a structural zero rather than having it floored at eps: Sinkhorn scaling
% needs full support, and the Bethe gap is a state-dependent lower bound that
% does not cancel when the estimates are normalized against each other.
%
% Reference:
%   H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical Monographs
%   14, Mathematical Association of America, 1963.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M, R] = size(L);
n = n(:)';
N = N(:)';

if nargin < 4, infset = []; end
if nargin < 5, lGN = []; end
if nargin < 6 || isempty(engine), engine = 'exact'; end

if numel(n) ~= M
    line_error(mfilename,'The occupancy vector has %d entries but L has %d rows.', numel(n), M);
end
if numel(N) ~= R
    line_error(mfilename,'The population vector has %d entries but L has %d columns.', numel(N), R);
end
if any(n < 0)
    line_error(mfilename,'The occupancy vector has a negative entry.');
end
if any(infset < 1 | infset > M)
    line_error(mfilename,'INFSET indexes a station outside 1..%d.', M);
end

% Infeasible occupancies are not an error: the caller sweeps a lattice.
if sum(n) ~= sum(N)
    pjoint = 0;
    lpjoint = -Inf;
    return
end

if isempty(lGN)
    isinfrow = false(1, M);
    isinfrow(infset) = true;
    Lq = L(~isinfrow, :);
    if any(isinfrow)
        Z = sum(L(isinfrow, :), 1);
    else
        Z = zeros(1, R);
    end
    [~, lGN] = pfqn_ca(Lq, N, Z);
end

if sum(N) == 0
    lpjoint = -lGN;
    pjoint = exp(lpjoint);
    return
end

% The expanded matrix is square of order SUM(N). 'spm' works on the
% unexpanded form, and building this would throw away the very property that
% makes it independent of the population.
if strcmpi(engine, 'spm')
    A = zeros(0, 0);
else
    A = sub_replicate(L, N, n);
end

if ~strcmpi(engine, 'exact')
    [i0, r0] = sub_firstzero(L, N, n);
    if ~isempty(i0)
        line_error(mfilename, ['The ''%s'' permanent engine cannot be applied: the demand of class %d ' ...
            'at station %d is zero, so the replicated matrix has no full support. ' ...
            'Use engine ''exact''.'], engine, r0, i0);
    end
end

switch lower(engine)
    case 'exact'
        F = pfqn_perm(A);
    case 'spm'
        % Never the expanded A: the saddle point is asymptotic in the column
        % multiplicities, which are the class populations themselves.
        [Ar, mr] = sub_replicate_rows(L, N, n);
        F = perm_spm(Ar, mr);
    case 'bethe'
        F = perm_bethe(A);
    case 'heur'
        F = perm_heur(A);
    case 'huberlaw'
        F = perm_huberlaw(A);
    case 'adapart'
        F = perm_adapart(A);
    otherwise
        line_error(mfilename,'Unrecognized permanent engine ''%s''. Use exact, spm, bethe, heur, huberlaw or adapart.', engine);
end

if F <= 0
    pjoint = 0;
    lpjoint = -Inf;
    return
end

lpjoint = log(F) - sum(factln(N)) - sum(factln(n(infset))) - lGN;
pjoint = exp(lpjoint);
end

function A = sub_replicate(L, N, n)
% Column R of L repeated N(R) times, then row I of that repeated n(I) times.
% A station holding no jobs and a class holding no jobs each drop out here,
% which is what makes a zero entry of the occupancy vector free of any special
% case: the result stays square of order sum(N).
M = size(L, 1);
Ak = [];
for r = 1:size(L, 2)
    if N(r) > 0
        Ak = [Ak, repmat(L(:,r), 1, N(r))]; %#ok<AGROW>
    end
end
A = zeros(0, size(Ak, 2));
for i = 1:M
    if n(i) > 0
        A((end+1):(end+n(i)), :) = repmat(Ak(i,:), n(i), 1);
    end
end
end

function [Ar, m] = sub_replicate_rows(L, N, n)
% Row I of L repeated n(I) times, keeping only the classes that hold jobs, with
% those classes' populations as the column multiplicities. This is the same
% matrix SUB_REPLICATE expands, one step earlier: PERM(AR, M) == PERM(A), and
% PERM_SPM wants the unexpanded form because its expansion is asymptotic in M.
% A class with no jobs is dropped rather than passed with multiplicity zero, so
% a zero demand in such a column cannot trip the full-support check.
keepc = find(N > 0);
m = N(keepc);
Lk = L(:, keepc);
Ar = zeros(0, numel(keepc));
for i = 1:size(L, 1)
    if n(i) > 0
        Ar((end+1):(end+n(i)), :) = repmat(Lk(i,:), n(i), 1);
    end
end
end

function [i0, r0] = sub_firstzero(L, N, n)
% First (station, class) whose zero demand actually reaches the replicated
% matrix. A class with no jobs or a station with no jobs contributes nothing,
% so its zeros are irrelevant.
i0 = [];
r0 = [];
for i = 1:size(L, 1)
    if n(i) == 0, continue; end
    for r = 1:size(L, 2)
        if N(r) > 0 && L(i,r) <= 0
            i0 = i;
            r0 = r;
            return
        end
    end
end
end
