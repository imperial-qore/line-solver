%{
%{
 % @file pfqn_lap.m
 % @brief Laplace approximation for normalizing constant.
%}
%}

%{
%{
 % @brief Laplace approximation for normalizing constant.
 % @fn pfqn_lap(L, N, Z)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @return logI Logarithm of normalizing constant approximation.
%}
%}
function [logI] = pfqn_lap(L,N,Z)
% [LOGI] = PFQN_LAP(L,N,Z)

Ntot = sum(N);
% f = @(x) 1-x + sum((N+x*L)./(Z+x*L));
% expv0 = fzero(f,1);
% logIv = log(Ntot) - sum(factln(N)); % coeff
% logIv = logIv -  + sum(N.*log(Z+L.*expv0)); % f(u0)
% logIv = logIv + 0.5*log(2*pi) - 0.5*log(sum((N./Ntot)./(Z./(Ntot.*L)+u0).^2)) ; % f''(u0)
% logIv = logIv - 0.5*log(Ntot);

Ntot = sum(N);
f = @(x) 1-sum(N.*L./(Z+Ntot*L.*x));
u0 = fzero(f,1,optimset('Display','off'));
if u0<0
    logI = NaN;
    return
end
%    f2 = sum((N./Ntot)./(Z./(Ntot.*L)+u0).^2);
%    logI = log(1-normcdf(0,u0,f2));
%    return
%end
if ~isfinite(u0)
    initSign = f(0.001)/abs(f(0.001));
    for u0 = 1e-4:1e-4:10
        fu0 = f(u0);
        if fu0/abs(fu0) ~= initSign
            break;
        end
    end
end
if u0<0
    logI = NaN;
    return
end
logI = log(Ntot) - sum(factln(N)); % coeff
logI = logI - Ntot*u0 + sum(N.*log(Z+L.*u0*Ntot)); % f(u0)
logI = logI + 0.5*log(2*pi) - 0.5*log(sum((N./Ntot)./(Z./(Ntot.*L)+u0).^2)) ; % f''(u0)
logI = logI - 0.5*log(Ntot);

% absf2 = sum((N./Ntot)./(Z./(Ntot.*L)+u0).^2); % |f''(u0)|
% f2 =  sum((N./Ntot)./(Z./(Ntot.*L)+u0).^2);
% f3 = -2 * sum((N./Ntot)./(Z./(Ntot.*L)+u0).^3);
% f4 = 6 * sum((N./Ntot)./(Z./(Ntot.*L)+u0).^4);
%
% logI3 = log(Ntot) - sum(factln(N)); % coeff
% logI3 = logI3 - Ntot*u0 + sum(N.*log(Z+L.*u0*Ntot)); % f(u0)
% logI3 = logI3 + (1/2)*log(2*pi/absf2) - (3/2)*log(Ntot);
% logI3 = logI3 + log(f4/(8*f2^2) - 5*f3^2/(24*f2^3));
end
