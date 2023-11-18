function key = hf1(marking,shift)
    m = length(marking);
    key = int32(0);
    slide = 0;
    for i=1:m
        key = bitxor(key, bitror(fi(marking(i)),slide).int32);
        slide = mod((slide + shift),32);
    end
end