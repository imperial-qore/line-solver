function [enabled_m, enabled_t] = enabledTransitions(F, INH, IM, P, M0)
[num_m, num_t] = size(F);
enabled_t = [];
enabled_m = [];
for t=1:num_t
    for m=1:num_m
        fi = cell2mat(F(m,t));
        if ~isempty(fi)
            in = cell2mat(INH(m,t));
            enabled = all(M0 >= fi & ~any(in>0 & in <= M0));
            if enabled
                if isempty(enabled_m)
                    enabled_m(end+1) = m;
                    enabled_t(end+1) = t;
                else
                    em = enabled_m(end);
                    et = enabled_t(end);
                    isImmed_1 = cell2mat(IM(m,t));
                    isImmed_2 = cell2mat(IM(em,et));
                    if isImmed_1 && isImmed_2
                        % Both mode is immediate
                        priority_1 = cell2mat(P(m,t));
                        priority_2 = cell2mat(P(em,et));
                        if priority_1 > priority_2
                            enabled_m = m;
                            enabled_t = t;
                        elseif priority_1 == priority_2
                            enabled_m(end+1) = m;
                            enabled_t(end+1) = t;
                        end
                    elseif isImmed_1 && ~isImmed_2
                        enabled_m = m;
                        enabled_t = t;
                    elseif ~isImmed_1 && ~isImmed_2
                        enabled_m(end+1) = m;
                        enabled_t(end+1) = t;
                    end
                end
            end
        end
    end
end
end