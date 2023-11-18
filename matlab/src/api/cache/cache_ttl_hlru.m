function prob=cache_ttl_hlru(lambda,m)
lambda = lambda(1,:,2:end);
lambda = reshape(lambda,size(lambda,2),size(lambda,3));
[n,h]=size(lambda); % n totalitems, h lists

%t = che_approx_T(lambda,m);  % FPI for the characteristic time of each list, a total of h lists
t = cache_t_hlru(lambda,m);
probh = zeros(n,1);   % steady state at list h
prob0 = zeros(n,1);   % steady state at list 0
temp1 = ones(n,1);
temp2 = zeros(n,1);

for k = 1:n
    for s = 1:h
        temp1(k) = temp1(k)*(1-exp(-lambda(k,s)*t(s)));
    end

    middtemp = 1;
    middtemp2 = 0;
    for l = 1:h-1
        for s = 1:l
            middtemp = middtemp*(1-exp(-lambda(k,s)*t(s)));
        end
        middtemp2 = middtemp2+middtemp;
    end
    temp2(k) = exp(-lambda(k,h)*t(h))*(1+middtemp2);

    probh(k) = temp1(k)/(temp1(k)+temp2(k));
    prob0(k) = 1/(temp1(k)+temp2(k))/exp(lambda(k,h)*t(h));
end
prob = [prob0 probh];
end