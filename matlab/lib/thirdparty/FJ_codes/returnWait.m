function [ wait_alpha, wait_Smat, prob_wait, alfa ]=returnWait( En1,pi0,T,phi,sum_Ajump )
%
% The function calculates the phase-type representation of the stationary
% waiting time of a K-node FJ queue
%

%% Section 6
alfa = -pi0/T; 
ds = length(phi); 

rhos = diag(phi*alfa)'/(alfa*phi);
En0 = alfa*sum_Ajump/sum(pi0);

prob_wait = (En0-1)/(En0-1+En1);
wait_alpha = prob_wait*rhos;

wait_Smat = zeros(ds,ds);
for i = 1:ds
    for j = 1:ds
        wait_Smat(i,j) = alfa(j)*T(j,i)/alfa(i); % Gij
    end
end

