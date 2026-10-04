function [jumps, n_dep] = ode_jumps_new(M, K, enabled, q_indices, P, Kic)
% [JUMPS, N_DEP] = ODE_JUMPS_NEW(M, K, MATCH, Q_INDICES, P, KIC)
%
% N_DEP is the number of leading columns that are service COMPLETIONS; the
% columns after them are intra-PH phase changes. ODE_RATE_BASE emits its rates
% in the same order, so a caller can tell which (station,class) an event is a
% completion of.

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
% Everything emitted so far is a service completion; what follows is an
% intra-PH phase change.
n_dep = size(jumps, 2);
for i = 1 : M   %state changes: "next service phase" transition
    for c = 1:K
        if enabled(i,c)
            xic = q_indices(i,c);
            % every source phase, the last included: bounding ki at Kic-1 is
            % valid only for an acyclic PH and drops the last row of D0 for a
            % general MAP or MMPP2, whose D0 is cyclic
            for ki = 1 : Kic(i,c)
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
