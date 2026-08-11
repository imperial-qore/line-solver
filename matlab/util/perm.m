function val = perm(A, m)
% Computes the permanent of matrix A applying computational savings if 
% some rows or columns are repeated

if nargin==1
    % Find unique columns and their indices
    [uniqueCols, ~, ic] = unique(A', 'rows', 'stable');

    % Count occurrences of each unique column
    m = histcounts(ic, 1:(size(uniqueCols,1)+1));

    % If no repeated columns (m is all ones), check for repeated rows
    if all(m == 1)
        [uniqueRows, ~, ir] = unique(A, 'rows', 'stable');
        % Count occurrences of each unique row
        m_rows = histcounts(ir, 1:(size(uniqueRows,1)+1));

        % If there are repeated rows, transpose and use row multiplicities
        if any(m_rows > 1)
            A = uniqueRows';
            m = m_rows;
        else
            % No repeated columns or rows, use original matrix
            A = uniqueCols';
        end
    else
        % There are repeated columns, use unique columns
        A = uniqueCols';
    end
end

R = length(m);
n = sum(m);
val = 0;
f = pprod(m);

while f >= 0
    term = (-1)^sum(f);

    % Multinomial coefficients
    for j=1:R
        term = term * nchoosek(m(j), f(j));
    end

    % Product term
    for i=1:n
        sumterm = 0;
        for k=1:R
            sumterm = sumterm + f(k) * A(i,k);
        end
        term = term * sumterm;
    end

    val = val + term;
    f = pprod(f, m);
end

val = (-1)^n * val;
end
