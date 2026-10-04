function prot = lqn_overtake_markov(phasetab, prvisit, srvresid, ycalls)
% PROT = LQN_OVERTAKE_MARKOV(PHASETAB, PRVISIT, SRVRESID, YCALLS)
%
% Overtaking probability for a phased server from the jump chain of
% Franks (1999), Ch. 5. An arrival that reaches the server while it is
% still executing phase 2 of an earlier request from the same client is
% an overtaking arrival; PROT is the probability of that event.
%
% The chain state is (p,x,k): p is the client phase, x the server phase,
% and k = 0,1,2 marks client execution, client dispatch, and server
% execution respectively (Franks 1999, Fig. 5.10). The chain has a
% product form and is solved column by column:
%
%   Eq. (3.1)  a phase of the client is cut into NSLICES = 1 + sum_j y_ijp
%              slices of host execution, each of mean SERVICE/NSLICES
%   Eq. (5.5)  Pr(OT_p | S_r,x), absorption into overtaking
%   Eq. (5.8)  Pr(NEXT_p | S_r,x), same expression with a_p replacing c_p
%   Eq. (5.6)  Pr{OT(x)}, the sum over start states of the column
%   Eq. (5.7)  Pr{S_p,2,0}, start probabilities of the first column
%
% This routine solves the first column (x = 2) for the case in which the
% calling entry is also the conditioning entry, so the start probabilities
% of Eq. (5.7) reduce to the phase carrying the call and Eq. (5.9) is not
% needed. Client phases beyond 2 and distinct conditioning entries are not
% built by SolverLN, so the columns they would generate are absent.
%
% Inputs
%   PHASETAB : (maxphase+1) x 5, row p+1 holds the parameters of client
%              phase p as [NSLICES SERVICE Y_IJ Y_IK T_K], where SERVICE is
%              the host residence of the phase, Y_IJ the calls to the
%              server task, Y_IK the calls to every other task, and T_K the
%              mean delay incurred at those other tasks. Row 1 (p = 0) is
%              the client think slice, so NSLICES = 1 and SERVICE = Z.
%   PRVISIT  : probability that the client visits this entry, 1 for a
%              reference task or a sole entry.
%   SRVRESID : residence time of the server phase under test, s_jx.
%   YCALLS   : (maxphase+1) vector, YCALLS(1) total calls to the server
%              task, YCALLS(p+1) the calls issued from client phase p.
%
% Output
%   PROT     : probability that an arrival finds the server in phase x.
%
% Reference: G. Franks, "Performance Analysis of Distributed Server
% Systems", PhD thesis, Carleton University, 1999, Sec. 5.4; published as
% G. Franks and M. Woodside, "Effectiveness of early replies in
% client-server systems", Perform. Eval. 36 (1999) 165-183.

maxphase = size(phasetab,1) - 1;
nphases = maxphase + 1;

% Per-phase coefficients a_p, b_p, c_p, d_p of Eq. (5.5) and Eq. (5.8)
ap = zeros(1,nphases); bp = zeros(1,nphases);
cp = zeros(1,nphases); dp = zeros(1,nphases);
for idx = 1:nphases
    p = idx - 1;
    if p == maxphase, pra = prvisit; else, pra = 1.0; end
    [ap(idx),bp(idx),cp(idx),dp(idx)] = local_jump_rates(srvresid, pra, ...
        phasetab(idx,1), phasetab(idx,2), phasetab(idx,3), ...
        phasetab(idx,4), phasetab(idx,5));
end

% Absorption probabilities of Eq. (5.5) and Eq. (5.8), for every start
% state (p,x,0) of the column and every conditioning phase r
prot_pr = zeros(nphases, nphases);
prnext_pr = zeros(nphases, nphases);
for p = 0:maxphase
    tail = 1.0;
    for r = 0:maxphase
        if r == p, continue; end
        tail = tail * local_ratio_b(bp(r+1), dp(r+1));
    end
    weight = 1.0 / local_absorb_denom(bp(p+1), dp(p+1), tail);
    r = p;
    while true
        prot_pr(p+1, r+1) = cp(p+1) * weight;    % Eq. (5.5)
        prnext_pr(p+1, r+1) = ap(p+1) * weight;  % Eq. (5.8)
        if r == 0, r = maxphase; else, r = r - 1; end
        weight = weight * local_ratio_b(bp(r+1), dp(r+1));
        if r == p, break; end
    end
end

% Eq. (5.7) with the calling entry as the conditioning entry: the start
% probability sits entirely on the phase that carries the call
prstart = zeros(1, maxphase);
if maxphase >= 1, prstart(maxphase) = 1.0; end

% Eq. (5.6) summed over start states, each weighted by the share of the
% flow to the server task that phase p contributes
prot = 0.0;
for p = 1:maxphase
    if phasetab(p+1,3) == 0, continue; end
    acc = 0.0;
    for r = 1:maxphase
        acc = acc + prstart(r) * prot_pr(p+1, r+1);
    end
    prot = prot + acc * (ycalls(1) / ycalls(p+1));
end
end

% -------------------------------------------------------------------------
function [ap,bp,cp,dp] = local_jump_rates(srvresid, pra, nslices, service, y_ij, y_ik, t_k)
% Coefficients of Eq. (5.5) and Eq. (5.8) from the seven jump probabilities
% q_{p,x,k} of Franks (1999), p. 107, written with mean times rather than
% rates, so mu_jx = 1/srvresid, mu_ip = 1/slice and mu_kp = 1/t_k.
ytot = y_ij + y_ik + 1.0;                         % Y_ip + 1
if nslices ~= 0, slice = service / nslices; else, slice = 0.0; end   % Eq. (3.1)

% q0, q3: race between the server phase and a reply from some other task
den = srvresid + t_k;
if ~isfinite(den)
    q0 = 0.0; q3 = 1.0;
elseif den ~= 0
    q0 = t_k / den; q3 = srvresid / den;
else
    q0 = 1.0; q3 = 0.0;
end

% q1, q5: race between the server phase and the client execution slice
den = srvresid + slice;
if ~isfinite(den)
    q1 = 0.0; q5 = 1.0;
elseif den ~= 0
    q1 = srvresid / den; q5 = slice / den;
else
    q1 = 1.0; q5 = 0.0;
end

q2 = y_ik / ytot;   % client calls some task other than the server
q4 = y_ij / ytot;   % client calls the server and overtakes
q6 = pra / ytot;    % client completes phase p and moves to phase p+1

ap = q5 + q1*q2*q0;
bp = q1*q6;
cp = q1*q4;
dp = q1*q2*q3;
end

function v = local_ratio_b(bp, dp)
% The factor b_y/(1 - d_y) of the products in Eq. (5.5) and Eq. (5.8)
v = bp / (1.0 - dp);
end

function v = local_absorb_denom(bp, dp, tail)
% The shared denominator of Eq. (5.5) and Eq. (5.8)
v = 1.0 - (bp * tail + dp);
end
