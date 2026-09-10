function n = ansim_roundpreserve(x, total)
% ANSIM_ROUNDPRESERVE  Largest-remainder rounding of the nonnegative vector X to
% integers summing exactly to TOTAL. Used to discretize the frozen field of the
% ansim surrogate when it is solved as an exact CTMC.
x = max(x(:).', 0);
n = zeros(size(x));
if total <= 0
    return;
end
if sum(x) <= 0
    n(1) = total;
    return;
end
x = x * total / sum(x);
n = floor(x);
short = total - sum(n);
if short > 0
    [~, ord] = sort(x - n, 'descend');
    n(ord(1:short)) = n(ord(1:short)) + 1;
end
end
