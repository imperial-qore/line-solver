function value = l1norm_cell(vectors)
value = 0.0;
for k = 1:numel(vectors)
    value = value + norm(vectors{k}, 1);
end
end
