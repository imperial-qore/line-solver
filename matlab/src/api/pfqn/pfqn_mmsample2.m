function [G,lG,lGr] = pfqn_mmsample2(L,N,Z,samples)
% [G,LOGG,LOGGR] = PFQN_MMSAMPLE2(L,N,Z,SAMPLES)

% Monte carlo sampling for normalizing constant of a repairmen model
% based on McKenna-Mitra integral form
R = length(N);
% Scale so that all coefficients are >=1.
scaleFactor = 1e-7 + min([L(:);Z(:)]); L = L/scaleFactor; Z = Z/scaleFactor;
nnzeros = 1:R;
c = 0.5;
v = [rand(1,ceil(c*samples)),logspace(0,5,ceil(samples*(1-c)))]; % sample more below the mean of the exponential
du = [0,diff(v)]';
u  = repmat(v',1,length(nnzeros));
ZL = log(repmat(Z+L(1,nnzeros),size(u,1),1).*u);
lG = du + -v' + ZL*N';
den = sum(factln(N));
lG = max(lG) - den; % get largest
lG = lG + sum(N)*log(scaleFactor); % rescale
for r=1:R
    lGr = lG - ZL(:,r);
    lGr(r) = max(lGr) - den + log(N(r)) + (sum(N)-1)*log(scaleFactor);
end
G = exp(lG);
end
