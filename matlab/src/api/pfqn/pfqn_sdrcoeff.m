%{
%{
 % @file pfqn_sdrcoeff.m
 % @brief Derived coefficients of a Krzesinski state-dependent routing structure.
%}
%}

%{
%{
 % @brief Derived coefficients of a Krzesinski state-dependent routing structure.
 % @fn pfqn_sdrcoeff(sdr)
 % @param sdr State-dependent routing structure.
 % @return c Derived coefficient structure.
%}
%}
function c = pfqn_sdrcoeff(sdr)
% C = PFQN_SDRCOEFF(SDR)
%
% Validates a state-dependent routing (SDR) structure and returns the derived
% coefficients of Krzesinski (1987), "Multiclass Queueing Networks with
% State-Dependent Routing", Performance Evaluation 7:125-143, eqs. (11)-(14).
%
% SDR describes the partition of the network into a subnetwork Q(V,V) subject
% to SDR and its complement Q(N-V,M-V), with Q(V,V) split into branches
% arranged in a hierarchy of nested subnetworks V_1 > V_2 > ... > V_T. Its
% fields are, following the paper's own indexing in which branch 1 is the
% complement M-V and the SDR branches are numbered 2..B:
%
%   sdr.entry       station index of the entry center e of Q(V,V)
%   sdr.departure   station index of the departure center d of Q(V,V)
%   sdr.branch      1xB cell, branch{b} = station indices of branch b (b>=2);
%                   branch{1} is unused and holds the complement implicitly
%   sdr.entryOf     1xB, entryOf(b) = station index of the entry center e(b)
%   sdr.departureOf 1xB, departureOf(b) = station index of d(b)
%   sdr.level       1xB, level(b) = the unique t with B_b in V_t - V_{t+1}
%   sdr.C           1xT, the coefficients C_t of eq. (11)
%   sdr.d           TxB, d(t,b) defined for 2<=b<=B and 1<=t<=level(b)
%
% The returned C has fields
%   c.T, c.B, c.level, c.C, c.d as above;
%   c.inA{t}        branch indices b with level(b) >= t, i.e. the set A_t of
%                   branches contained in the subnetwork Q(V_t,V_t)
%   c.Dtt(t)        D_tt = sum_{b in A_t} d(t,b), eq. (14)
%   c.Dprev(t)      D_{t-1,t} = sum_{b in A_t} d(t-1,b), eq. (14), t>1
%   c.mmax(b)       largest m_b with delta_{level(b),b}(m_b) >= 0
%   c.vmax(t)       largest v_t with omega_tt(v_t) >= 0 and, for t>1, with
%                   omega_{t-1,t}(v_t) >= 0
%
% The bounds c.mmax and c.vmax are the population constraints that SDR imposes
% on branches and on subnetworks (Sec. 2.5); they are consequences of the
% routing coefficients, not independent inputs. A bound is Inf when the
% corresponding C_t is nonnegative, in which case SDR imposes no constraint.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isstruct(sdr)
    line_error(mfilename,'The SDR structure must be a struct.');
end
for f = {'entry','departure','branch','entryOf','departureOf','level','C','d'}
    if ~isfield(sdr,f{1})
        line_error(mfilename,sprintf('The SDR structure is missing field ''%s''.',f{1}));
    end
end

B = numel(sdr.branch);
T = numel(sdr.C);
if B < 2
    line_error(mfilename,'An SDR structure must declare at least one branch (branch indices start at 2).');
end
if T < 1
    line_error(mfilename,'An SDR structure must declare at least one level of subnetwork nesting.');
end
level = sdr.level(:)';
if numel(level) ~= B
    line_error(mfilename,'SDR level must have one entry per branch index.');
end
if any(level(2:end) < 1) || any(level(2:end) > T) || any(level(2:end) ~= round(level(2:end)))
    line_error(mfilename,sprintf('SDR branch levels must be integers in 1..%d.',T));
end
if size(sdr.d,1) < T || size(sdr.d,2) < B
    line_error(mfilename,sprintf('SDR coefficient matrix d must be at least %dx%d.',T,B));
end

% The nesting V_1 > V_2 > ... > V_T must be strict: every level must carry at
% least one branch, otherwise two hierarchically adjacent subnetworks coincide
% and the ratio omega_{t-1,t}/omega_tt in eq. (10) is not the paper's.
for t = 1:T
    if ~any(level(2:end) == t)
        line_error(mfilename,sprintf('SDR level %d carries no branch: the subnetwork nesting must be strict.',t));
    end
end

% Branches must be disjoint, must not contain the entry or departure center of
% Q(V,V), and must contain their own declared entry and departure centers.
seen = [];
for b = 2:B
    sb = sdr.branch{b}(:)';
    if isempty(sb)
        line_error(mfilename,sprintf('SDR branch %d is empty.',b));
    end
    if any(ismember(sb,seen))
        line_error(mfilename,'SDR branches must be mutually disjoint.');
    end
    if any(sb == sdr.entry) || any(sb == sdr.departure)
        line_error(mfilename,'The entry and departure centers of Q(V,V) must not belong to any branch.');
    end
    if ~any(sb == sdr.entryOf(b)) || ~any(sb == sdr.departureOf(b))
        line_error(mfilename,sprintf('The entry and departure centers of SDR branch %d must belong to that branch.',b));
    end
    seen = [seen, sb]; %#ok<AGROW>
end

inA = cell(1,T);
Dtt = zeros(1,T);
Dprev = zeros(1,T);
for t = 1:T
    inA{t} = find(level >= t & (1:B) >= 2);
    Dtt(t) = sum(sdr.d(t,inA{t}));
    if t > 1
        Dprev(t) = sum(sdr.d(t-1,inA{t}));
    end
end

% Population bounds implied by the requirement that delta and omega stay
% nonnegative (Sec. 2.5). With C_t < 0 these are finite and are exactly the
% concurrency constraints SDR places on branches and subnetworks.
mmax = Inf(1,B);
for b = 2:B
    t = level(b);
    if sdr.C(t) < 0
        mmax(b) = floor(sdr.d(t,b) / (-sdr.C(t)));
    end
end
vmax = Inf(1,T);
for t = 1:T
    if sdr.C(t) < 0
        vmax(t) = floor(Dtt(t) / (-sdr.C(t)));
    end
    if t > 1 && sdr.C(t-1) < 0
        vmax(t) = min(vmax(t), floor(Dprev(t) / (-sdr.C(t-1))));
    end
end

c = struct();
c.T = T;
c.B = B;
c.level = level;
c.C = sdr.C(:)';
c.d = sdr.d;
c.inA = inA;
c.Dtt = Dtt;
c.Dprev = Dprev;
c.mmax = mmax;
c.vmax = vmax;
c.entry = sdr.entry;
c.departure = sdr.departure;
c.branch = sdr.branch;
c.entryOf = sdr.entryOf;
c.departureOf = sdr.departureOf;
end
