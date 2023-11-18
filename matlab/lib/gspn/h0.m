function key = h0(marking, N)
%     key = double(mod(mod(hf1(marking,17),37),N)+1);
    key = double(mod(mod(hf3(marking),37),N)+1);
end