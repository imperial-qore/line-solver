%{ @file cache_ttl_lrua.m
 %  @brief Computes steady-state probabilities for TTL-LRU cache with arrivals
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes steady-state probabilities for TTL-LRU cache with arrivals
 %
 % @details
 % This function computes the steady-state probability distribution for a
 % TTL-LRU cache system with multiple users, items, and cache levels.
 %
 % @par Syntax:
 % @code
 % prob = cache_ttl_lrua(lambda, R, m)
 % prob = cache_ttl_lrua(lambda, R, m, seed)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rates per user per item per list
 % <tr><td>R<td>Routing probability structure
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>seed<td>(Optional) Random seed for initialization
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>prob<td>Steady-state probability distribution
 % </table>
%}
function prob=cache_ttl_lrua(lambda, R, m, seed)
if nargin<4
    seed = 23000;
end
rng(seed,'twister');

u=size(lambda,1); % number of users
n=size(lambda,2); % number of items
h=size(lambda,3)-1; % number of lists

fun = @ttl_tree_time;

% random seed to generate the initial value for x
range_left = 0; range_right = 10;
x = (range_right-range_left).*rand(1,h) + range_left;
options = optimoptions('fsolve','MaxIter',1e5,'MaxFunEvals',1e6,'Display','off');
[listtime,~,~,~] = fsolve(fun,x,options);
[~,ssprob] = ttl_tree_time(listtime);
prob = ssprob;

    function [F,ssprob,capadiff] = ttl_tree_time(x)
        steadystateprob = zeros(n,h+1);
        randprob = zeros(n,h+1);
        avgtime = zeros(n,h+1);
        cdiff = zeros(1,h);
        capa = zeros(1,h);
        rpdenominator = zeros(1,n);
        % the probability of each item at each list
        for i = 1:n % for all items
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
                randprob(i,dtchain(a)) = steadystateprob(i,dtchain(a))*avgtime(i,dtchain(a))/rpdenominator(i); % random time , prob of item i in list j
            end
        end
        ssprob = randprob;

        for l = 1:h
            capa(l) = sum(randprob(:,l+1));
            cdiff(l) = m(l)-capa(l);
        end
        capadiff = cdiff;
        F = cdiff;
    end

end
