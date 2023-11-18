function tangible = isTangible(F, INH, IM, M0)
[num_m, num_t] = size(F);
for t=1:num_t
    for m=1:num_m
        fi = F{m,t};
        if ~isempty(fi)
            in = INH{m,t};
            enabled = all(M0 >= fi & ~any(in>0 & in <= M0));
            if enabled
                isImmed = IM{m,t};
                if isImmed
                    tangible = false;
                    return
                end
            end
        end
    end
end
tangible = true;
end