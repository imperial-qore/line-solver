function dx = ode_pnorm(x, Phi, Mu, PH, M, K, enabled, q_indices, rt, Kic, nservers, w, sched_id, pstar)
% DX = ODE_PNORM(X, PHI, MU, PH, M, K, ENABLED, Q_INDICES, RT, KIC, NSERVERS, W, SCHED_ID, PSTAR)
%
% ODE derivative function using p-norm smoothing for processor-share constraint.
% Based on Ruuskanen et al., PEVA 151 (2021).
%
% This provides a smoother approximation than softmin, improving ODE stability
% for stiff problems.
%
% @param pstar Smoothing parameter vector (one per station) or scalar
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

dx = 0*x;

% If pstar is scalar, expand to per-station
if isscalar(pstar)
    pstar = pstar * ones(M, 1);
end

for i = 1:M
    p = pstar(i); % p-norm parameter for this station

    switch sched_id(i)
        case SchedStrategy.INF
            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p);
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
                                            if j~=i
                                                rate = Phi{i}{c}(kic) * Mu{i}{c}(kic) * rt((i-1)*K+c,(j-1)*K+l) * pie(kjl);
                                                dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1) * rate;
                                                dx(xjl+kjl-1) = dx(xjl+kjl-1) + x(xic+kic-1) * rate;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

        case SchedStrategy.EXT
            % Source node with external arrivals - p-norm smoothing not needed
            % (arrivals are independent of queue state)
            % Phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p);
                                dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate;
                                dx(xic+kic_p-1) = dx(xic+kic_p-1) + x(xic+kic-1)*rate;
                            end
                        end
                    end
                end
            end
            % Service completions (departures from source)
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
                                            if j~=i
                                                rate = Phi{i}{c}(kic) * Mu{i}{c}(kic) * rt((i-1)*K+c,(j-1)*K+l) * pie(kjl);
                                                dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1) * rate;
                                                dx(xjl+kjl-1) = dx(xjl+kjl-1) + x(xic+kic-1) * rate;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

        case SchedStrategy.PS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );

            % Compute p-norm smooth factor: ghat = 1 / (1 + (ni/c)^p)^(1/p)
            if ni > 0 && nservers(i) > 0
                ghat = pnorm_smooth(ni, nservers(i), p);
            else
                ghat = 1;
            end

            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p) * ghat;
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
                                            rate = Phi{i}{c}(kic) * Mu{i}{c}(kic) * rt((i-1)*K+c,(j-1)*K+l) * pie(kjl) * ghat;
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

        case SchedStrategy.FCFS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );
            wni = GlobalConstants.FineTol;
            for c=1:K
                for kic = 1 : Kic(i,c)
                    if enabled(i,c)
                        xic = q_indices(i,c);
                        w(c,kic) = -1/PH{i}{c}{1}(kic,kic);
                        wni = wni + w(c,kic)*x(xic + kic - 1);
                    end
                end
            end

            % Compute p-norm smooth factor
            if ni > 0 && nservers(i) > 0
                ghat = pnorm_smooth(ni, nservers(i), p);
            else
                ghat = 1;
            end

            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p) * ghat * w(c,kic) / wni;
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
                                            rate = rate * ghat * w(c,kic) / wni;
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

        case SchedStrategy.DPS
            w(i,:) = w(i,:)/sum(w(i,:));
            % Compute weighted queue length for rate normalization
            wni = GlobalConstants.FineTol;
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                wni = wni + w(i,k) * sum(x(idxIni:idxEnd));
            end
            % Compute total queue length for ghat
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );

            % Compute p-norm smooth factor using total queue length
            if ni > 0 && nservers(i) > 0
                ghat = pnorm_smooth(ni, nservers(i), p);
            else
                ghat = 1;
            end

            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p) * ghat * w(i,c) / wni;
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
                                            rate = rate * ghat * w(i,c) / wni;
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
