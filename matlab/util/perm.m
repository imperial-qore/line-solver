function val = perm(A, m)
% VAL = PERM(A)      permanent of the square matrix A
% VAL = PERM(A, M)   permanent of the matrix whose column J is column J of A
%                    repeated M(J) times, so that SUM(M) == SIZE(A,1)
%
% Computational savings are applied when rows or columns are repeated.

% Conditioning is a floating-point notion, so symbolic input keeps the caller's
% orientation and its supplied multiplicities unchanged.
numericInput = isnumeric(A);

if nargin==2 && numericInput
    % Materialize the expansion so both entry points share the orientation
    % choice below; without this the two-argument form is pinned to whichever
    % orientation the caller happened to build.
    A = repelem(A, 1, m(:)');
end

if numericInput
    % Grouping repeated rows means expanding the transpose, which leaves the
    % permanent unchanged but can raise the largest inclusion-exclusion term by
    % many orders of magnitude. Orientation is chosen first, grouping second;
    % see _kb/03-api-layer.md.
    if perm_conditioning(A') < perm_conditioning(A)
        A = A';
    end

    % Find unique columns of the chosen orientation and their multiplicities
    [uniqueCols, ~, ic] = unique(A', 'rows', 'stable');
    m = histcounts(ic, 1:(size(uniqueCols,1)+1));
    A = uniqueCols';
elseif nargin==1
    % Symbolic single-argument call: every column is its own group
    m = ones(1, size(A,2));
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
