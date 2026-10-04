function [ D0, D1, p ] = map2( H2lambda1, H2lambda2, H2p, decay)

%Heindl: "INVERSE CHARACTERIZATION OF HYPEREXPONENTIAL MAP(2)S"

% approximate the seperate arrival rates and probability
% r1, r2, r3: first three moments
r1 = H2p/H2lambda1 + (1-H2p)/H2lambda2;
r2 = (2*H2p/(H2lambda1^2) + 2*(1-H2p)/(H2lambda2^2))/(factorial(2));
r3 = (6*H2p/(H2lambda1^3) + 6*(1-H2p)/(H2lambda2^3))/(factorial(3));

% the first moment should be same for all three distributions
h2 = (r2-r1^2)/(r1^2);
h3 = (r1*r3-r2^2)/(r1^4);
a = h2;
d = h3 + h2^2;
b = d - a;
c = b^2 + 4*a^3;

if decay >= 0
    lambda1 = 2*h2/(r1*(2*a+b-sqrt(c)));
    lambda2 = 2*h2/(r1*(2*a+b+sqrt(c)));
    p = 0.5 + b/(2*sqrt(c));
    D0 = [-lambda1,0; 0,-lambda2];
    D1 = (decay-1)*[-lambda1; -lambda2]*[p,1-p] - decay*D0;
elseif decay < 0
    x = (b+sqrt(c))/(2*a^2);
    y = (-b+sqrt(c))/(2*a^2);
    gamma1 = -x*(a^2*y+d)/(a*(y+1));
    gamma2 = y*(a^2*x-d)/(a*(x-1));
    if b >= 0
        gamma3 = -y/x;
        
    elseif b < 0
        gamma3 = -x/y;
        
    end
    lower_bound = max(max(gamma1,gamma2),gamma3);
    decay = ceil(lower_bound*10)/10;
    fprintf('The smallest decay should be %0.2f, now decay is %0.2f \n', lower_bound, decay);
    
    D0 = [-(2*a+b+sqrt(c)),0; 0,-(2*a+b-sqrt(c))]/(2*r1*h3);
    D11 = (2*a+b+sqrt(c))*(1+b/sqrt(c)+decay*(1-b/sqrt(c)));
    D12 = (2*a+b+sqrt(c))*(1-b/sqrt(c))*(1-decay);
    D21 = (2*a+b-sqrt(c))*(1+b/sqrt(c))*(1-decay);
    D22 = (2*a+b-sqrt(c))*(1-b/sqrt(c)+decay*(1+b/sqrt(c)));
    D1 = [D11, D12; D21, D22]/(4*r1*h3);
    p = x/(x+y);
end
end