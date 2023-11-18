function D=qmle(Q,N,Z)
[M,R] = size(Q);
D = 0*Q;
for i=1:M
    for j=1:R
        D(i,j)= Q(i,j) / (N(j)-sum(Q(:,j)))* Z(j) / (1+ sum(Q(i,:))- Q(i,j)/N(j));
    end
end
end