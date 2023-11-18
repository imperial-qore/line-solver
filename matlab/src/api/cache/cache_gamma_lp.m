function [gamma,u,n,h]=cache_gamma_lp(lambda,R)
u=size(lambda,1); % number of users
n=size(lambda,2); % number of items
h=size(lambda,3)-1; % number of lists

gamma=zeros(n,h);
for i = 1:n % for all items
    for j = 1:h % for all levels
        % compute gamma(i,j)
        Rvi = 0*R{1,i};
        for v=1:u
            Rvi = Rvi + R{v,i};
        end
        Pij =[1+j];
        pr_j = par(Rvi, 1+j);
        while ~isempty(pr_j)
            Pij =[pr_j, Pij];
            pr_j = par(Rvi, pr_j);
        end        
        if isempty(Pij)   % debug
        %if length(Pij) < 2
            gamma(i,j)=0;
        else
            gamma(i,j)=1;
            for li = 2:length(Pij) % for all levels up to the current one
                y = 0;
                l_1 = Pij(li-1);
                l = Pij(li);
                for v = 1:u % for all streams
                    for t=1:l_1
                        y=y + lambda(v, i, t) * R{v,i}(t, l);
                    end
                end
                gamma(i,j)=gamma(i,j)*y;
            end
        end
    end
end

end

function parent = par(R, j)
% finds the parent of j according to the access probabilities in R
    parent = find(R(1:(j-1),j));
    if length(parent) > 1
        line_error(mfilename,'A cache has a list with more than one parent, but the structure must be a tree.');
    end
end