% Parameters
n = 7;  
m = [1,1,3];  
h = length(m);
lambda = [49, 49, 49, 49, 7, 1, 1]/205;

% Initial condition (randomly initialized probabilities)
rng(23000, 'twister');
x0 = rand(n, 1+ h);
for i=1:n
    x0(i,:) = zeros(1,1+h);
    x0(i,1) = 1;
end    

% Time span
tspan = [0 10000];

% Solve ODE
[t, x] = ode23s(@(t, x) cache_rrm_meanfield_ode(t, x, lambda, m, n, h), tspan, x0);

x=reshape(x(end,:),[n,1+h])
missrate  = lambda*x(:,1)
missratio = (lambda*x(:,1))/sum(lambda)

