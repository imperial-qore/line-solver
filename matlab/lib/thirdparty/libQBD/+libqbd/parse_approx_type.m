function max_degree = parse_approx_type(approx_type)
raw_value = uint32(approx_type);
max_degree = double(bitand(raw_value, uint32(2^24 - 1)));
max_degree = min(max_degree, libqbd.get_max_factor());
end
