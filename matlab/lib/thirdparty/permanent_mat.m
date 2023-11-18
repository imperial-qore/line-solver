function p = permanent_mat( A )
% Computes the permament of square matrix, A
%   the entries are assumed to be real or complex

% Written by Brian K. Butler, 2009-2016

[m, n]=size(A);

if (m > n)
    error('Matrix must be have # rows <= # columns.  Error inside permament_mat() %d x %d\n',m,n) 
end

if (m == 0),
    p = 1;
elseif (m == 1)
    p = sum(A);
elseif (m == n) % square
    p = loc_permanent_mat_sq( A , m );
else   % non-square
    p = loc_permanent_mat_rect( A , m , n );
    
end

return

% local function or 'subfunction' for SQUARE MATRIX
function p = loc_permanent_mat_sq( A , m )

if (m == 2),
    p = A(1,1)*A(2,2) + A(1,2)*A(2,1);
else
    p = 0;
    for col = 1:m,
        %compind = [1:(ind-1) (ind+1):m];  % complementary indices up to m
        if A(1,col)~=0, %<-- allows complex
            p = p + A(1,col) * loc_permanent_mat_sq(A(2:m,[1:(col-1) (col+1):m]),m-1);
        end
    end
end

% local function or 'subfunction' for RECTANGULAR MATRIX
function p = loc_permanent_mat_rect( A , m , n )

if (m == 2),
    % from Henryk Minc book "Permanents" 1978, Addison-Wesley, pp 2-3.
    p = prod(sum(A,2)) - sum(prod(A,1));  
else
    p = 0;
    for col = 1:n,
        %compind = [1:(ind-1) (ind+1):n];  % complementary indices up to n
        if A(1,col)~=0, %<-- allows complex
            p = p + A(1,col) * loc_permanent_mat_rect(A(2:m,[1:(col-1) (col+1):n]), m-1, n-1);
        end
    end
end