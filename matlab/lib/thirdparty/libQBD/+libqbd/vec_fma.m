function a = vec_fma(a, b, cons)
n = min(numel(a), numel(b));
for k = 1:n
    a{k} = a{k} + cons * b{k};
end

for k = (n + 1):numel(b)
    a{end + 1} = cons * b{k};
end
end
