function e = enorm(matrix)
    % Computes the euclidean norm of a matrix
    e = 0;
    for i=1:size(matrix, 1)
        for j=1:size(matrix, 2)
            e = e + matrix(i, j)^2;
        end
    end
    e = sqrt(e);
end