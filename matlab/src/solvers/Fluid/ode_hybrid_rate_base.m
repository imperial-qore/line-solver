function [rateBase, eventIdx] = ode_hybrid_rate_base(Phi, Mu, PH, M, K, enabled, q_indices, rt, Kic, sched_id, all_jumps)
% [RATEBASE, EVENTIDX] = GETRATEBASE(PHI, MU, PH, M, K, MATCH, Q_INDICES, P, KIC, SCHED_ID, ALL_JUMPS)

% Phi{i}{c}(ki): probability of service completion in phase ki of station i
%               for class c jobs

rateBase = zeros(size(all_jumps,2),1);
eventIdx = zeros(size(all_jumps,2),1);
rateIdx = 0;

for i = 1 : M   %state changes from departures in service phases 2...
    for c = 1:K
        if enabled(i,c)
            for j = 1 : M
                for l = 1:K
                    pie = map_pie(PH{j}{l});
                    if rt((i-1)*K+c,(j-1)*K+l) > 0
                        for kic = 1 : Kic(i,c)
                            for kjl = 1 : Kic(j,l)
                                rateIdx = rateIdx + 1;
                                switch sched_id(i)
                                    case SchedStrategy.ID_FCFS
                                        rateBase(rateIdx) = 0;
                                    otherwise
                                        rateBase(rateIdx) = Phi{i}{c}(kic) * Mu{i}{c}(kic) * rt((i-1)*K+c,(j-1)*K+l) * pie(kjl);
                                end
                                eventIdx(rateIdx) = q_indices(i,c) + kic - 1;
                            end
                        end
                    end
                end
            end
        end
    end
end

for i = 1 : M   %state changes from "next service phase" transition in phases 2...
    for c = 1:K
        if enabled(i,c)
            for kic = 1 : (Kic(i,c) - 1)
                for kicp = 1 : Kic(i,c) % (ki+1), PH distribution
                    if kicp ~= kic
                        rateIdx = rateIdx + 1;
                        switch sched_id(i)
                            case SchedStrategy.ID_FCFS
                                rateBase(rateIdx) = 0;
                            otherwise
                                
                                rateBase(rateIdx) = PH{i}{c}{1}(kic,kicp);
                        end
                        eventIdx(rateIdx) = q_indices(i,c) + kic - 1;
                    end
                end
            end
        end
    end
end

end % getRateBase

