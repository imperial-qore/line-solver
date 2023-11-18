function dx = ode_hybrid(x, Phi, Mu, PH, M, K, enabled, q_indices, rt, Kic, nservers, w, sched_id, all_jumps, rateBase, eventIdx)
% RATES = ODE_HYBRID(X, M, K, Q_INDICES, KIC, NSERVERS, W, STRATEGY, RATEBASE, EVENTIDX)

dxrates = x; % basic vector valid for INF and PS case min(ni,nservers(i))=ni
for i = 1:M
    switch sched_id(i) % source
        case SchedStrategy.ID_INF
            % do nothing
        case SchedStrategy.ID_EXT  %EXT
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                dxrates(idxIni) = 1-sum(x(idxIni+1:idxEnd)); % not needed for idxIni+1:idxEnd as rates is initialized equal to x
            end
        case SchedStrategy.ID_PS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );
            if ni > nservers(i) % case  min = ni handled by rates = x
                dxrates(idxIni:idxEnd) = x(idxIni:idxEnd)/ni * nservers(i);
            end
        case SchedStrategy.ID_FCFS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            dxrates(idxIni:idxEnd) = 0;
        case SchedStrategy.ID_DPS %DPS
            w(i,:) = w(i,:)/sum(w(i,:));
            ni = mean(w(i,:));
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                ni = ni + sum( w(i,k)*x(idxIni:idxEnd) );
            end
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                dxrates(idxIni:idxEnd) = w(i,k)*x(idxIni:idxEnd)/ni * nservers(i); % not needed for idxIni+1:idxEnd as rates is initiliazed equal to x
            end
    end
end
dxrates = dxrates(eventIdx);
dxrates = rateBase.*dxrates;
dx = all_jumps * dxrates;

% obtain dx for state-dependent nodes
for i = 1:M
    switch sched_id(i)
        case SchedStrategy.ID_FCFS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            dx(idxIni:idxEnd) = 0;
            ni = sum( x(idxIni:idxEnd) );
            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic, kic_p);
                                if ni > nservers(i)
                                    rate = rate * nservers(i) /ni;
                                end
                                dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate;
                                dx(xic+kic_p-1) = dx(xic+kic_p-1) + x(xic+kic-1)*rate;
                            end
                        end
                    end
                end
            end
            % service completions
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for j = 1 : M
                        for l = 1:K
                            xjl = q_indices(j,l);
                            if enabled(j,l)
                                pie = map_pie(PH{j}{l});
                                if rt((i-1)*K+c,(j-1)*K+l) > 0
                                    for kic = 1 : Kic(i,c)
                                        for kjl = 1 : Kic(j,l)
                                            rate = Phi{i}{c}(kic) * Mu{i}{c}(kic) * rt((i-1)*K+c,(j-1)*K+l) * pie(kjl);
                                            if ni > nservers(i)
                                                rate = rate * nservers(i) /ni;
                                            end
                                            dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate;
                                            dx(xjl+kjl-1) = dx(xjl+kjl-1) + x(xic+kic-1)*rate;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
    end
end
end