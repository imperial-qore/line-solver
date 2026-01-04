function [W,rhohat]=qsys_gigk_approx_kingman(lambda,mu,ca,cs,k)
% W=QSYS_GIG1_UBND_KINGMAN(LAMBDA,MU,CA,CS,K)
W = (ca^2+cs^2)/2 * (qsys_mmk(lambda,mu,k)-1/mu)+1/mu;
rhohat = W*lambda/(1+W*lambda); % so that M/M/1 formulas still hold
end
