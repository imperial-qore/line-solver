%{
%{
 % @file pfqn_mu_ms.m
 % @brief Compute load-dependent rates for m identical c-server FCFS stations.
%}
%}

%{
%{
 % @brief Compute load-dependent rates for m identical c-server FCFS stations.
 % @fn pfqn_mu_ms(N, m, c)
 % @param N Maximum population.
 % @param m Number of identical stations.
 % @param c Number of servers per station.
 % @return mu Load-dependent service rate vector (1xN).
%}
%}
function mu = pfqn_mu_ms(N,m,c)
% calculates the load-dependent rate of m identical c-server FCFS stations
mu = zeros(1,N);
g = zeros(m,N+1); %table
for n=0:N
    for i=1:m
        g(i,1+n) = gnaux(n,i,c,g);
    end
end
for n=1:N
    mu(n) = g(m,1+(n-1))/g(m,1+n); % 1+ is to handle indexing from 1 instead than 0
end
end

function gn = gnaux(n,m,c,g)
if n==0
    gn = 1;
else
    if m==1
        gn = 1/prod(min(1:n,c*ones(1,n)));
    else
        gn = 0;
        for k=0:n
            a = min(1:k,c*ones(1,k)); 
            b = 1/g(m-1,1+n-k); 
            gn = gn + 1/(prod(a)*prod(b)); % prod([])=1
        end
    end
end
end