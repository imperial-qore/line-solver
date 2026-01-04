function W=polling_qsys_gated(arvMAPs,svcMAPs,switchMAPs)
% W=polling_qsys_gated(arvMAPs,svcMAPs,switchMAPs)
%
% Exact mean waiting time solution of a polling system with open arrivals.
% All queues use gated service.
%
% Takagi, ACM Computing Surveys, Vol. 20, No. 1, March 1988, eq (20)
%
% Example:
% W=polling_qsys_gated({map_exponential(1/0.6),map_exponential(1/0.2)},{map_exponential(1),map_exponential(1)},{map_exponential(1),map_exponential(1)})

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
r=sum(r1);
% station time method Ferguson and Aminetzah 1985 as reported by Takagi
W=gated(n,rho,delta2,lambda,b2,rho1,r);
end

function W=gated(n,rho,delta2,lambda,b2,rho1,r)
lst1 = [];
lst2 = [];

for i = 1:n
    for j = 1:n
        if i > j
            t1 = zeros(1,n^2);
            for m = i:n
                t1((j-1)*n + m) = -1;
            end
            for m = 1:j-1
                t1((j-1)*n + m) = -1;
            end
            for m = j:i-1
                t1((m-1)*n + j) = -1;
            end
            t1((i-1)*n + j) = 1 / rho1(i);
            lst1 = [lst1; t1];
            lst2 = [lst2; 0];
        elseif j > i
            t1 = zeros(1,n^2);
            for m = i:j-1
                t1((j-1)*n + m) = -1;
            end
            for m = j:n
                t1((m-1)*n + j) = -1;
            end
            for m = 1:i-1
                t1((m-1)*n + j) = -1;
            end
            t1((i-1)*n + j) = 1 / rho1(i);
            lst1 = [lst1; t1];
            lst2 = [lst2; 0];
        else
            t1 = zeros(1,n^2);
            t1((i-1)*n + i) = t1((i-1)*n + i) + 1;
            for m = 1:n
                if i ~= m
                    t1((i-1)*n + m) = -rho1(i);
                end
            end
            for m = 1:n
                t1((m-1)*n + i) = t1((m-1)*n + i)-(rho1(i)^2);
            end
            lst1 = [lst1; t1];
            temp = delta2(i) + lambda(i) * b2(i) * r / (1 - rho);
            lst2 = [lst2; temp];
        end
    end
end
final = linsolve(lst1, lst2);
W = zeros(1, n);

for i = 1:n
    temp = (1 + rho1(i)) * r / (2 * (1 - rho));
    sum = 0;
    for j = 1:n
        if i ~= j
            sum = sum + final((i-1)*n + j);
        end
    end
    sum = sum * (1 / rho1(i));
    for j = 1:n
        sum = sum + final((j-1)*n + i);
    end
    temp = temp + (1 - rho) * (1 + rho1(i)) * sum / (2 * r);
    %disp(temp);
    W(i) = temp;
end
end
