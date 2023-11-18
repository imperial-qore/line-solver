function COV = mtrace_cov(T, A)

C = max(A);
N = length(A);

COV = cell(C, C);

for c1 = 1:C
    for c2 = 1:C
        X0c1v = zeros(N-1,1);
        X1c2v = zeros(N-1,1);
        for i = 1:(N-1)
            X0c1 = 0;
            X1c2 = 0;
            if A(i) == c1
                X0c1 = T(i);
            end
            if A(i+1) == c2
                X1c2 = T(i+1);
            end
            X0c1v(i) = X0c1;
            X1c2v(i) = X1c2;
        end
        COV{c1,c2} = cov(X0c1v, X1c2v);
    end
end

end