function Y = swapcell(P)
% Given a cell array of identically sized matrices, e.g.
% P = cell(2,2);
% P{1,1} = [0,1,0; 0,0,1; 0,0,0];
% P{1,2} = [0,0,0; 0,0,0; 1,0,0];
% P{2,1} = [0,0,0; 0,0,0; 1,0,0];
% P{2,2} = [0,1,0; 0,0,1; 0,0,0];
% this method swaps station and class indexes, returning a cell array
% where Y{i,j}(r,s) = P{r,s}(i,j)

M = length(P{1,1});
K = length(P);
Y= cellzeros(M,M,K,K);
for i=1:M
    for j=1:M
        for r=1:K
            for s=1:K
                Y{i,j}(r,s) = P{r,s}(i,j);
            end
        end
    end
end
end