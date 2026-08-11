function value = function_of_dist(dist, func)
value = 0.0;
for k = 1:numel(dist)
    weights = func(k - 1, numel(dist{k}));
    value = value + sum(libqbd.as_row_vector(dist{k}) .* libqbd.as_row_vector(weights));
end
end
