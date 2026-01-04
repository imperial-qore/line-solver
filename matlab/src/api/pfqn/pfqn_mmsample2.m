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
v = [rand(1,ceil(c*samples)),logspace(0,5,ceil(samples*(1-c)))]; % sample more below the mean of the exponential
du = [0,diff(v)]';
u  = repmat(v',1,R);
ZL = log(repmat(Z+L(1,1:R),size(u,1),1).*u);
lG = du + -v' + ZL*N';
lG = max(lG) - sum(factln(N)); % get largest
lG = lG + sum(N)*log(scaleFactor); % rescale
G = exp(lG);
end
