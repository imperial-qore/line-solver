function W=polling_qsys_1limited(arvMAPs,svcMAPs,switchMAPs)
% W=polling_qsys_1limited(arvMAPs,svcMAPs,switchMAPs)
%
% Exact mean waiting time solution of a polling system with open arrivals.
% All queues use 1-limited service.
%
% O. J. Boxma and B. Meister. Waiting-time approximations for
% cyclic-service systems with switch-over times. SIGMETRICS
% /PERFORMANCE '86, page 254–262, New York, NY, USA, 1986.
%
% Example:
% W=polling_qsys_1limited({map_exponential(1/0.6),map_exponential(1/0.2)},{map_exponential(1),map_exponential(1)},{map_exponential(1),map_exponential(1)})

n = length(arvMAPs); % number of classes
for i=1:n
    lambda(i) = map_lambda(arvMAPs{i});
    b(i) = map_mean(svcMAPs{i});
    b2(i) = map_moment(svcMAPs{i},2);
    rho1(i) = lambda(i)*b(i);
    r1(i) = map_mean(switchMAPs{i});
    delta2(i) = map_var(switchMAPs{i});
end
rho=sum(rho1);
R=sum(r1);
% station time method Ferguson and Aminetzah 1985 as reported by Takagi
W=onelimited(n,rho,delta2,lambda,b2,rho1,R);
end

function W=onelimited(n,rho,delta2,lambda,b2,rho1,R)
W = zeros(1,n);
for i=1:n
    W(i)=(1-rho+rho1(i))/(1-rho-lambda(i)*R);
    W(i)=W(i)*(1-rho)/((1-rho)*rho+sum(rho1.^2));
    W(i)=W(i)*(rho/(2*(1-rho))*sum(lambda*b2(:))+rho*sum(delta2)/2/R + R/(2*(1-rho))*(rho1(i)*(1+rho1(i))));
    end
end
