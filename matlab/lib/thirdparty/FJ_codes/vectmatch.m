function i = vectmatch(row,matrix)

% find the corresponding row in the matrix
% here we assume there is one and only one 'row' in matrix

[m,n] = size(matrix);
% walk from end of matrix array and search for row starting with row.

for outer = 1:m
    for inner = 1:n+1
        if inner == n+1
            i = outer;
            return;
        else
            if matrix(outer,inner) ~= row(inner)
                break;
            end
        end
    end
end

