function [gamma,u,n,h]=cache_gamma(lambda,R)
u=size(lambda,1); % number of users
n=size(lambda,2); % number of items
h=size(lambda,3)-1; % number of lists

gamma=zeros(n,h);
for i = 1:n % for all items
    for j = 1:h % for all levels
        %%
        % compute gamma(i,j)        
        
        G = digraph(R{1,i}); % we take the topology from stream 1 as must be identical across streams?
        Pj = G.shortestpath(1,j);
        
        gamma(i,j)=sum(lambda(:,i,1+0)); % level 1        
        for li = 2:length(Pj) % for all levels up to the current one
            y = 0;        
            l_1 = Pj(li-1);
            l = Pj(li);
            for v = 1:u % for all streams
                y=y + lambda(v, i, 1 + l_1) * R{v,i}(l_1, l);                
            end
            gamma(i,j)=gamma(i,j)*y;
        end
    end
end
end