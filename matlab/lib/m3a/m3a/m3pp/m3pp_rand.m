function MMAP=m3pp_rand(order,classes)
MAP = mmpp_rand(order);
MMAP{1} = MAP{1};
MMAP{2} = MAP{2};
for i=1:order
    p = rand(1,classes);
    p = p / sum(p);
    for j=1:order
        for c=1:classes
            MMAP{2+c} = MMAP{2}*p(c);
        end
    end
end
end
