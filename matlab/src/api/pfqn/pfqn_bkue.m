%{
%{
 % @file pfqn_bkue.m
 % @brief Birman-Kogan uniform (van der Waerden) expansion for a single chain.
%}
%}

%{
%{
 % @brief Birman-Kogan uniform (van der Waerden) expansion for a single chain.
 %
 % Birman and Kogan (Stochastic Models 8(3):543-563, 1992), Section 4. The
 % plain saddle point loses accuracy once the saddle approaches the dominant
 % pole of the integrand, which is the regime where the station holding that
 % pole saturates. The van der Waerden uniform expansion keeps the pole and
 % the saddle in one formula through the complementary error function, and so
 % stays accurate on both sides of the crossing.
 %
 % The published formula carries two typographical defects that the paper's
 % own Table 3 settles; see the notes in the code below.
 %
 % @fn pfqn_bkue(L, N, Z)
 % @param L Service demand vector (stations x 1), single class.
 % @param N Population (scalar).
 % @param Z Think time (default: 0).
 % @return G Normalizing constant.
 % @return lG Logarithm of the normalizing constant.
%}
%}
function [G,lG] = pfqn_bkue(L,N,Z)
% [G,LG] = PFQN_BKUE(L,N,Z)
%
% Uniform expansion of the single chain normalizing constant.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(Z)
    Z = 0;
end
Z = sum(Z(:));
L = L(:);
L(L<=0) = [];
if isempty(N) || N == 0
    G = 1; lG = 0;
    return
end
if isempty(L)
    lG = N*log(Z) - gammaln(N+1); % the delay alone
    G = exp(lG);
    return
end
% The dominant pole is the slowest station, and it is a pole rather than part
% of the exponent only when it stands alone against a group: a station with an
% identical twin is one of the paper's groups and belongs in the exponent, and
% without any group at all there is no M_j >> 1 against which a lone station
% is an O(1) factor, so the expansion degenerates to the plain saddle point.
[dmax,ipole] = max(L);
tolL = GlobalConstants.FineTol*max(1,dmax);
hasGroup = numel(L) > 1 && any(abs(diff(sort(L))) <= tolL);
hasPole = nnz(abs(L-dmax) <= tolL) == 1 && hasGroup;
if hasPole
    D = L; D(ipole) = [];
    zp = 1/dmax;
else
    D = L;
    zp = inf;
end
z0 = local_saddle(D,N,Z);
h2 = local_d2(z0,D,N);
h3 = local_d3(z0,D,N);
t2 = (1/z0 + h3/(6*h2))/sqrt(2*pi*h2);
if ~isfinite(zp)
    % No pole to keep out of the exponent, so every station is in it. The
    % expansion degenerates to the plain saddle point, and the third
    % derivative term goes with the pole it corrects: it belongs to the
    % regular part of the integrand once the pole has been subtracted, and
    % keeping it with no pole to subtract misstates lG by nats.
    lG = local_h(z0,D,N,Z) - log(z0) - 0.5*log(2*pi*h2);
    G = exp(lG);
    return
end
b2 = local_h(zp,D,N,Z) - local_h(z0,D,N,Z); % the paper's M*b^2, always >= 0
if zp >= z0 % saddle before the pole
    lG = local_h(z0,D,N,Z) + log(0.5*erfcx(sqrt(b2)) + t2);
else % the pole has been crossed and its residue leads
    lG = local_h(zp,D,N,Z) + log(1 - 0.5*erfc(sqrt(b2)) + t2*exp(-b2));
end
G = exp(lG);
end

function f = local_h(z,D,N,Z)
f = Z*z - N*log(z) - sum(log(1 - D*z));
end

function g = local_d1(z,D,N,Z)
g = Z - N/z + sum(D./(1 - D*z));
end

function h = local_d2(z,D,N)
h = N/z^2 + sum(D.^2./(1 - D*z).^2);
end

function h = local_d3(z,D,N)
h = -2*N/z^3 + 2*sum(D.^3./(1 - D*z).^3);
end

function z = local_saddle(D,N,Z)
if isempty(D)
    z = N/Z;
    return
end
hi = 1/max(D);
z = 0.5*hi;
for it = 1:200
    g = local_d1(z,D,N,Z);
    if abs(g) <= 1e-14*max(1,N)
        break
    end
    dz = -g/local_d2(z,D,N);
    alpha = 1;
    while z + alpha*dz <= 0 || z + alpha*dz >= hi
        alpha = alpha/2;
        if alpha < 1e-14
            break
        end
    end
    if alpha < 1e-14
        break
    end
    z = z + alpha*dz;
end
end
