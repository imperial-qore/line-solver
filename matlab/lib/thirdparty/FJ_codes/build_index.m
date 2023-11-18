function indexes = build_index(m,cr)
total_dim = nchoosek(cr+m-1,cr);
indexes = zeros(total_dim,m);
indexes(1,1) = cr;
for row = 2 : total_dim
    k = find(indexes(row-1,:)>0, 1);
    if k < m
        indexes(row,:) = indexes(row-1,:);
        indexes(row,k+1) = indexes(row,k+1)+1;
        indexes(row,1) = indexes(row,k)-1;
        indexes(row,2:k) = 0;                     
    end
    
end


