%{
%{
 % @file pfqn_mmsample2.m
 % @brief Monte Carlo sampling for repairman models using McKenna-Mitra form.
%}
%}

%{
%{
 % @brief Monte Carlo sampling for repairman models using McKenna-Mitra form.
 % @fn pfqn_mmsample2(L, N, Z, samples)
 % @param L Service demand vector.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param samples Number of samples.
 % @return G Normalizing constant estimate.
 % @return lG Logarithm of normalizing constant.
%}
%}
function [G,lG] = pfqn_mmsample2(L,N,Z,samples)
% [G,LG] = PFQN_MMSAMPLE2(L,N,Z,SAMPLES)

% Monte carlo sampling for normalizing constant of a repairmen model
% based on McKenna-Mitra integral form
R = length(N);
% Scale so that all coefficients are >=1.
scaleFactor = 1e-7 + min([L(:);Z(:)]); 
L = L/scaleFactor; 
Z = Z/scaleFactor;
c = 0.5;
% The quadrature nodes must be SORTED: v mixes uniform draws with a logspace
% grid, so diff(v) on the raw concatenation is sign-indefinite and the panel
% widths below would be negative.
v = sort([rand(1,ceil(c*samples)),logspace(0,5,ceil(samples*(1-c)))]); % sample more below the mean of the exponential
du = [v(1),diff(v)]'; % panel widths, first panel covering [0,v(1)]
u  = repmat(v',1,R);
% McKenna-Mitra: G = 1/prod(N_r!) * int_0^inf exp(-u) prod_r (Z_r + L_r u)^N_r du.
% The integrand is (Z_r + L_r*u), NOT (Z_r + L_r)*u -- the latter is a different
% function that happens to have the same value at u = 1.
ZL = log(repmat(Z,size(u,1),1) + repmat(L(1,1:R),size(u,1),1).*u);
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
lterms = log(du) - v' + ZL*N';
lmax = max(lterms);
lG = lmax + log(sum(exp(lterms-lmax))) - sum(factln(N));
lG = lG + sum(N)*log(scaleFactor); % rescale
G = exp(lG);
end
