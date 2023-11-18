function prob=cache_ttl_tree(lambda, R, m)

global ssprob capadiff

u=size(lambda,1); % number of users
n=size(lambda,2); % number of items
h=size(lambda,3)-1; % number of lists

fun = @time;
rng(0,'twister');
range_left = 0;
range_right = 10;
x = (range_right-range_left).*rand(1,h) + range_left;    % random seed to generate the initial value for x
% x = ones(1,h)*0.5;
%options = optimoptions('fsolve','Display','iter','MaxIter',1e8,'MaxFunEvals',1e8,'StepTolerance',1e-6,'FunctionTolerance',1e-6);
%[options.OptimalityTolerance,options.FunctionTolerance,options.StepTolerance]
%options = optimoptions('fsolve','Display','iter','MaxIter',1e6,'MaxFunEvals',1e6);
% options = optimset('PlotFcns',@optimplotfval);
% [listtime,fval,exitflag,output] = fminsearch(fun,x,options);
%[listtime,fval,exitflag,output] = fsolve(fun,x,options);
[listtime,fval,exitflag,output] = fsolve(fun,x);
% sym('time_%d', [1 h]);
% listtime = vpasolve(fun);
%capadiff;
prob = ssprob;

function F = time(x)
steadystateprob = zeros(n,h+1);  
randprob = zeros(n,h+1);
avgtime = zeros(n,h+1);
cdiff = zeros(1,h);
capa = zeros(1,h);
rpdenominator = zeros(1,n);
% the probability of each item at each list
for i = 1:n % for all items
    %sym('time_%d', [1 4]);
    transmatrix = zeros(h+1, h+1);
    for j = 1:h+1
        leafnode = find(R{1,i}(j,:));
        for k = leafnode
            if j == 1
                transmatrix(j,k) = R{1,i}(j,k);
            else                
                transmatrix(j,k) = (1-exp(-lambda(1,i,j)*x(j-1)))*R{1,i}(j,k);
            end
            if j ~= k 
                transmatrix(k,j) = exp(-lambda(1,i,k)*x(k-1));      
            end           
        end       
    end
    missconnection = find(all(transmatrix==0));
    dtchain = setdiff(1:h+1, missconnection);
    transmatrix(missconnection,:)=[];
    transmatrix(:,missconnection)=[];   % remove the unused nodes in the transfer matrix
    dtmcprob = dtmc_solve(transmatrix);   % solution of dtmc, i.e., prob of item i in list j, 1*h
    for a = 1:numel(dtchain)
        steadystateprob(i,dtchain(a)) = dtmcprob(a);
        if dtchain(a)>1
            avgtime(i,dtchain(a)) = (1-exp(-lambda(1,i,dtchain(a))*x(dtchain(a)-1)))/lambda(1,i,dtchain(a));% average time of item i spent in l , l>=1
        else
            avgtime(i,dtchain(a)) = 1/lambda(1,i,dtchain(a));
        end
        rpdenominator(i) = rpdenominator(i)+steadystateprob(i,dtchain(a))*avgtime(i,dtchain(a));% denominator for the probability of item i at node j at a random time
    end
    for a = 1:numel(dtchain)
        %if dtchain(a)>1
            randprob(i,dtchain(a)) = steadystateprob(i,dtchain(a))*avgtime(i,dtchain(a))/rpdenominator(i); % random time , prob of item i in list j
        %end
    end
end
% x
% steadystateprob
% for q = 1:n
%     randprob(q,1) = 1-sum(randprob(q,2:end));
% end
ssprob = randprob;

for l = 1:h
    capa(l) = sum(randprob(:,l+1));
    cdiff(l) = m(l)-capa(l);
end
capadiff = cdiff;
F = cdiff;
end

end

