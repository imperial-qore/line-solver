function pie = mmap_pie(MMAP)
% For each class c, computes the stationary probability of the DTMC
% embedded at the restart intants after an arrival of class c.
% Input:
% - MMAP: the MMAP of order n with m classes
% Output:
% - pie:  an mxn matrix, where the c-th row is the solution for class c

symbolic = map_issym(MMAP);

n = size(MMAP{1},1);
m = size(MMAP,2)-2;

if ~symbolic
    pie = zeros(m,n);
else
    pie = sym(zeros(m,n));
end

for c = 1:m
    Pc = (-MMAP{1}-MMAP{2}+MMAP{2+c}) \ MMAP{2+c};
    if ~symbolic
        A = Pc'-eye(n);
        A(end,:) = ones(1,n);
        b = zeros(n,1);
        b(n) = 1;
    else
        A = Pc'-sym(eye(n));
        A(end,:) = sym(ones(1,n));
        b = sym(zeros(n,1));
        b(n) = 1;
    end
    pie(c,:) = (A \ b)';
end

end