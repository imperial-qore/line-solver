function res = mam_bgchain_station(Da0, Da1, alpha_s, T, A, esup, nservers, gref, Kmax, options)
% RES = MAM_BGCHAIN_STATION(DA0, DA1, ALPHA_S, T, A, ESUP, NSERVERS, GREF, KMAX, OPTIONS)
%
% Solves the open classes of one station as a level-dependent QBD MODULATED by
% the background chain: level = number of open jobs held by the station, phase
% = (arrival MAP phase, environment state, service phase).
%
% The station is an MAP/PH/c queue whose server is shared with the closed jobs.
% With k open and e closed jobs present the open aggregate completes at rate
%
%   phi(k,e) = min(k+e, c) * k/(k+e)
%
% times the phase-type completion rate of one busy server, i.e. the open class
% receives the share k/(k+e) of the min(k+e,c) busy servers. The dependence on
% k is what makes the QBD level-dependent; the dependence on e is what makes it
% modulated. The c parallel servers are collapsed into a single phase-type
% process scaled by phi, which is exact for exponential service at any c and
% for phase-type service at c = 1, and approximates the multiset of in-service
% phases otherwise -- the same collapse SOLVER_MAM_LDQBD documents.
%
% The environment is level-dependent too, and for the same reason. A lumped
% transition that LOWERS the closed occupancy of this station is a closed
% completion here, so it carries the closed share min(e+k,c)*e/(e+k) of the
% server; the background chain built it at the averaged share GREF, and level k
% rescales it by the ratio of the two. A transition that RAISES the occupancy is
% an arrival from elsewhere and is left alone. Without this the closed jobs
% would drain at their mean-field rate however long the open queue is, and the
% positive correlation between the two occupancies -- the very thing a congested
% station produces -- would be lost.
%
% The level space is truncated at KMAX. An arrival at the top level is lost but
% still advances the arrival phase, so the arrival process keeps its exact
% marginal and autocorrelation and only the queue tail is cut; RES.ploss
% reports the probability mass sitting at the top level.
%
% Inputs
%   DA0, DA1  (ma x ma)  aggregate open arrival MAP at the station
%   ALPHA_S   (1 x ms)   initial vector of the open service phase-type law
%   T         (ms x ms)  subgenerator of the open service phase-type law
%   A         (me x me)  environment generator, either the LUMPED chain from
%                        MAM_BGCHAIN_ENV or the unlumped one from
%                        MAM_BGCHAIN_ENVFULL; this function only needs to know
%                        how many closed jobs each of its states holds here
%   ESUP      (1 x me)   closed jobs each environment state stands for. It is
%                        strictly increasing for the lumped chain and REPEATS
%                        for the unlumped one, where many chain states hold the
%                        same number of jobs at this station
%   NSERVERS  scalar     number of servers
%   GREF      (1 x me)   closed capacity share A was built at
%   KMAX      scalar     truncation level of the open queue
%
% Output (struct RES)
%   .QLen   mean number of open jobs at the station
%   .Util   mean fraction of the servers held by open jobs
%   .Tput   open departure rate
%   .ploss  stationary probability of the truncation level
%   .plev   (1 x KMAX+1) level distribution of the open queue
%   .penv   (1 x ne)  stationary probability of holding each DISTINCT number of
%           closed jobs; states that hold the same number are folded together,
%           which is the identity for the lumped environment
%   .cshare (1 x ne)  E[min(e+k,c)*e/(e+k) | e], the closed capacity share the
%           background chain reads back
%   .esup   (1 x ne)  the distinct entries of ESUP, ascending: the index of
%           .penv and .cshare
%
% See also SOLVER_MAM_BGCHAIN, MAM_BGCHAIN_ENV, LDQBD_R, LDQBD_PI.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ma = size(Da0, 1);
me = numel(esup);
ms = numel(alpha_s);
alpha_s = alpha_s(:)';
t = -T * ones(ms, 1);

Ima = eye(ma);
Ime = eye(me);
Ims = eye(ms);

Kmax = max(1, round(Kmax));

Q0 = cell(Kmax, 1);
Q1 = cell(Kmax + 1, 1);
Q2 = cell(Kmax, 1);

% Environment split by WHAT A TRANSITION DOES TO THE CLOSED OCCUPANCY HERE:
% a departure from this station (which the open level throttles), an arrival to
% it (which it does not), and -- when the environment is the unlumped background
% chain -- a move that leaves this station's occupancy alone, which is a closed
% job travelling between two OTHER stations and is likewise untouched. The
% split is by ESUP rather than by triangle because ESUP REPEATS in that case;
% with the lumped environment, whose occupancies are distinct and sorted, the
% masks reduce to tril/triu and Aeq is empty, so the two agree entry for entry.
[~, eorder] = sort(esup(:)');           % stable: ties keep the chain's own order
A     = A(eorder, eorder);
esup  = esup(eorder);
gref  = gref(eorder);
isdown = bsxfun(@gt, esup', esup);      % row s -> col s' lowers the occupancy
isup   = bsxfun(@lt, esup', esup);
Adown = A .* isdown;
Aup   = A .* isup;
Aeq   = A - Adown - Aup;                % same occupancy, diagonal included
Aeq(1:me+1:end) = 0;                    % rebuilt at every level

Q1{1} = kron(Da0, Ime) + kron(Ima, bgchain_env_level(Aeq, Aup, Adown, esup, gref, 0, nservers));
Q0{1} = kron(kron(Da1, Ime), alpha_s);

phiae = cell(Kmax, 1);
for k = 1:Kmax
    share = bgchain_share(k, esup, nservers);
    phiae{k} = repmat(share, 1, ma);        % ordering (a outer, e inner)
    Ak = bgchain_env_level(Aeq, Aup, Adown, esup, gref, k, nservers);
    Q1{k+1} = kron(kron(Da0, Ime), Ims) + kron(kron(Ima, Ak), Ims) ...
        + kron(diag(phiae{k}), T);
    if k < Kmax
        Q0{k+1} = kron(kron(Da1, Ime), Ims);
    end
    if k == 1
        Q2{1} = kron(diag(phiae{1}), t);
    else
        Q2{k} = kron(diag(phiae{k}), t * alpha_s);
    end
end
% truncation: an arrival at the top level is lost, its phase transition is kept
Q1{Kmax+1} = Q1{Kmax+1} + kron(kron(Da1, Ime), Ims);

ldopts = struct('epsilon', options.tol, 'maxIter', options.iter_max, 'verbose', false);
R = ldqbd_R(Q0, Q1, Q2, ldopts);
[plev, pcell] = ldqbd_pi(R, Q0, Q1, Q2, ldopts);

plev = max(plev, 0);
plev = plev / sum(plev);

QLen = (0:Kmax) * plev(:);
Util = 0;
Tput = 0;
penv = zeros(1, me);
gacc = zeros(1, me);        % sum_k P(level k, environment e) * closed share
for k = 0:Kmax
    pk = max(pcell{k+1}(:)', 0);
    if k == 0
        marg = sum(reshape(pk, me, ma), 2)';
    else
        marg = sum(reshape(sum(reshape(pk, ms, me * ma), 1), me, ma), 2)';
        Util = Util + pk * kron(phiae{k}(:), ones(ms, 1));
        Tput = Tput + pk * kron(phiae{k}(:), t);
    end
    penv = penv + marg;
    gacc = gacc + marg .* bgchain_closed_share_vec(esup, k * ones(1, me), nservers);
end
Util = Util / nservers;
penvsum = sum(penv);
if penvsum > 0
    penv = penv / penvsum;
    gacc = gacc / penvsum;
end
% The caller reads the share back per NUMBER of closed jobs held, so states that
% hold the same number are folded here. With the lumped environment there is
% nothing to fold and this is the identity.
[eu, ~, eidx] = unique(esup);
penvu = accumarray(eidx(:), penv(:), [numel(eu), 1])';
gaccu = accumarray(eidx(:), gacc(:), [numel(eu), 1])';
cshare = zeros(1, numel(eu));
nzenv = penvu > GlobalConstants.Zero;
cshare(nzenv) = gaccu(nzenv) ./ penvu(nzenv);

res = struct('QLen', QLen, 'Util', Util, 'Tput', Tput, 'ploss', plev(Kmax+1), ...
    'plev', plev, 'penv', penvu, 'cshare', cshare, 'esup', eu);
end

function v = bgchain_share(k, esup, nservers)
% Share of the server capacity that k open jobs hold when esup closed jobs are
% also present: min(k+e,c) busy servers times the open fraction k/(k+e).
tot = k + esup(:)';
v = zeros(1, numel(tot));
nz = tot > 0;
v(nz) = min(tot(nz), nservers) .* (k ./ tot(nz));
end

function v = bgchain_closed_share_vec(esup, kvec, nservers)
% Share of the server capacity that the esup closed jobs hold, each entry
% against its own open occupancy KVEC. The mirror image of BGCHAIN_SHARE, and
% the factor the lumped closed departure rate is proportional to.
esup = esup(:)';
tot = esup + kvec(:)';
v = zeros(1, numel(tot));
nz = tot > 0;
v(nz) = min(tot(nz), nservers) .* (esup(nz) ./ tot(nz));
end

function Ak = bgchain_env_level(Aeq, Aup, Adown, esup, gref, k, nservers)
% Environment generator seen at open level k: the closed departures FROM THIS
% STATION are rescaled from the mean-field share GREF to the share they hold
% against k open jobs; the closed arrivals to it and the moves that leave its
% occupancy alone are unchanged; the diagonal is rebuilt.
g = bgchain_closed_share_vec(esup, k * ones(1, numel(esup)), nservers);
ratio = ones(1, numel(g));
nz = gref > 0;
ratio(nz) = g(nz) ./ gref(nz);
Ak = Aeq + Aup + diag(ratio) * Adown;
Ak = Ak - diag(sum(Ak, 2));
end
