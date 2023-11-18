function jumps = ode_jumps_new(M, K, enabled, q_indices, P, Kic)
% JUMPS = ODE_JUMPS_NEW(M, K, MATCH, Q_INDICES, P, KIC, STRATEGY)

jumps = []; %returns state changes triggered by all the events
jump = zeros( sum(sum(Kic)), 1 );
for i = 1 : M   %state changes from departures in service phases 2...
    for c = 1:K
        if enabled(i,c)
            xic = q_indices(i,c); % index of  x_ic
            for j = 1 : M
                for l = 1:K
                    if P((i-1)*K+c,(j-1)*K+l) > 0
                        xjl = q_indices(j,l); % index of x_jl
                        for ki = 1 : Kic(i,c) % job can leave from any phase in i
                            for kj = 1 : Kic(j,l) % job can start from any phase in j
                                jump = 0*jump; % reuse same vector for efficiency
                                jump(xic+ki-1) = jump(xic+ki-1) - 1; %type c in stat i completes service
                                jump(xjl+kj-1) = jump(xjl+kj-1) + 1; %type c job starts in stat j
                                jumps = [jumps jump;];
                            end
                        end
                    end
                end
            end
        end
    end
end
for i = 1 : M   %state changes: "next service phase" transition
    for c = 1:K
        if enabled(i,c)
            xic = q_indices(i,c);
            for ki = 1 : (Kic(i,c) - 1)
                for kip = 1:Kic(i,c)
                    if ki~=kip
                        jump = 0*jump; % reuse same vector for efficiency
                        jump(xic+ki-1) = jump(xic+ki-1) - 1;
                        jump(xic+kip-1) = jump(xic+kip-1) + 1;
                        jumps = [jumps jump;];
                    end
                end
            end
        end
    end
end
end % ode_jumps_new()
