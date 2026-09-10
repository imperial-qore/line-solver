function value = get_one_div_by_factor(n_plus_1)
if n_plus_1 > (libqbd.get_max_factor() - 1)
    value = 0.0;
    return;
end

value = exp(-gammaln(double(n_plus_1) + 2.0));
end
