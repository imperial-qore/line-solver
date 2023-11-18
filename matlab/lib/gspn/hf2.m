function key = hf2(marking, shift1, shift2)
    m = length(marking);
    key = int32(0);
    slide1 = 0;
    slide2 = 0;
    sum = 0;
    for i=1:m
       sum = sum + marking(i);
       key = bitxor(key, bitror(fi(marking(i)),slide1).int32);
       key = bitxor(key, bitror(fi(sum),slide2).int32);
       slide1 = mod((slide1 + shift1),32);
       slide2 = mod((slide2 + shift2),32);
    end
end