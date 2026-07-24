function [ S, s, sigma ] = get_serviceDis( mean, choice, CX2 )

% return the PH representation (s,S) of the required distribution
% mean: mean rate of the distribution
% choice: the targeted distribution, where 1 = Exp, 2 = HE2, 3 = ER2

if (choice == 1)
    % Exponential distribution:Exp
    S = -mean;
    s = mean;
    sigma = 1;
    
elseif (choice == 2)
    % 2-phase Hyper-exponential distribution: HE2
    [ lambda1,lambda2,a ] = estimateH2( mean,CX2 );
    S = [-lambda1,0; 0,-lambda2];
    s = [lambda1; lambda2];
    sigma = [a,1-a];
    
elseif (choice == 3)
    % 2-phase Erlang distribution: ER2
    v = 2*mean;
    S = [-v,v;0,-v];
    s = [0;v];
    sigma = [1,0];
end



