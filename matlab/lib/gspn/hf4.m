function total = hf4(marking)
total = 0;
m = length(marking);
for i=1:m
    total = 2 * total;
    for j=1:m-i
       total = total + marking(j);
    end
end
if total < 0
    total = intmax + total;
end
end