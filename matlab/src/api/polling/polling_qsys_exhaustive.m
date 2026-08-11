function W=polling_qsys_exhaustive(arvMAPs,svcMAPs,switchMAPs)
% W=polling_qsys_exhaustive(arvMAPs,svcMAPs,switchMAPs)
%
% Exact mean waiting time solution of a polling system with open arrivals.
% All queues use exhaustive service.
%
% Takagi, ACM Computing Surveys, Vol. 20, No. 1, March 1988, eq (15)
%
% Example:
% W=polling_qsys_exhaustive({map_exponential(1/0.6),map_exponential(1/0.2)},{map_exponential(1),map_exponential(1)},{map_exponential(1),map_exponential(1)})


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
W=exhaustive(n,rho,delta2,lambda,b2,rho1,R);
end

function W=exhaustive(n,rho,delta2,lambda,b2,rho1,r)
lst1 = [];
lst2 = [];

for i = 1:n
    for j = 1:n
        if i > j
            t1 = zeros(1, n^2);
            for m = (i + 1):n
                t1((j - 1) * n + m) = t1((j - 1) * n + m) - 1;
            end
            for m = 1:(j - 1)
                t1((j - 1) * n + m) = t1((j - 1) * n + m) - 1;
            end
            for m = j:(i - 1)
                t1((m - 1) * n + j) = t1((m - 1) * n + j) - 1;
            end
            t1((i - 1) * n + j) = t1((i - 1) * n + j) + ((1 - rho1(i)) / rho1(i));
            lst1 = [lst1; t1];
            t2 = 0;
            lst2 = [lst2; t2];
        elseif j > i
            t1 = zeros(1, n^2);
            for m = (i + 1):(j - 1)
                t1((j - 1) * n + m) = t1((j - 1) * n + m) - 1;
            end
            for m = j:n
                t1((m - 1) * n + j) = t1((m - 1) * n + j) - 1;
            end
            for m = 1:(i - 1)
                t1((m - 1) * n + j) = t1((m - 1) * n + j) - 1;
            end
            t1((i - 1) * n + j) = t1((i - 1) * n + j) + ((1 - rho1(i)) / rho1(i));
            lst1 = [lst1; t1];
            t2 = 0;
            lst2 = [lst2; t2];
        else % case r_ii
            t1 = zeros(1, n^2);
            t1((i - 1) * n + i) = t1((i - 1) * n + i) + 1;
            for m = 1:n
                if i ~= m
                    t1((i - 1) * n + m) = t1((i - 1) * n + m) - (rho1(i) / (1 - rho1(i)));
                end
            end
            lst1 = [lst1; t1];
            temp = 0;
            if i>1
                temp = temp + ((delta2(i - 1)) / (1 - rho1(i))^2);
            else
                temp = temp + ((delta2(n)) / (1 - rho1(i))^2);
            end
            temp = temp + (lambda(i) * b2(i) * r * (1 - rho1(i)) / ((1 - rho) * ((1 - rho1(i))^3)));
            lst2 = [lst2; temp];
        end
    end
end

final = lst1\lst2;
%disp(final);

W = zeros(1, n);

for i = 1:n
    temp = 0;
    temp = temp + (lambda(i) * b2(i)) / (2 * (1 - rho1(i)));
    temp = temp + r * (1 - rho1(i)) / (2 * (1 - rho));
    sum = 0;
    for j = 1:n
        if i ~= j
            sum = sum + final((i - 1) * n + j);
        end
    end
    sum = sum * ((1 - rho1(i)) / rho1(i));
    if i>1
        sum = sum + delta2(i - 1);
    else
        sum = sum + delta2(n);
    end
    sum = sum / (r * (1 - rho1(i)) * 2 / (1 - rho));
    temp = temp + sum;
    %disp(temp);
    W(i) = temp;
end
end
