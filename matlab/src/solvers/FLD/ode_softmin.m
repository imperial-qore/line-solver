function dx = ode_softmin(x, Phi, Mu, PH, M, K, enabled, q_indices, rt, Kic, nservers, w, sched_id, alpha)
% RATES = ODE_SOFTMIN(X, M, K, Q_INDICES, KIC, NSERVERS, W, SCHED_ID, ALPHA)

% Variance of ODE_STATEDEP with softmin function replacing the min function

dx = 0*x;
for i = 1:M
    switch sched_id(i) % source
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
        case SchedStrategy.EXT  %EXT
            % For open models, use options.method='matrix' or options.method='pnorm' instead
            line_error(mfilename,'Softmin ODE method does not support open models. Use options.method=''matrix'' or ''pnorm'' instead.');
        case SchedStrategy.PS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );
            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p);
                                if ni > nservers(i)
                                    dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate* nservers(i) /ni;
                                    dx(xic+kic_p-1) = dx(xic+kic_p-1) + x(xic+kic-1)*rate* nservers(i) /ni;
                                else
                                    dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate;
                                    dx(xic+kic_p-1) = dx(xic+kic_p-1) + x(xic+kic-1)*rate;
                                end
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
                                                rate = 1/ni * nservers(i) * rate;
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
            
            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p);
                                rate = rate * softmin(ni,nservers(i),alpha) * w(c,kic) /wni;
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
                                            rate = rate * softmin(ni,nservers(i),alpha) * w(c,kic) /wni;
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
        case SchedStrategy.DPS %DPS
            w(i,:) = w(i,:)/sum(w(i,:));
            % Compute weighted queue length for rate normalization
            wni = GlobalConstants.FineTol;
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                wni = wni + w(i,k) * sum(x(idxIni:idxEnd));
            end
            % Compute total queue length for capacity comparison
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );
            % phase changes
            for c = 1:K
                if enabled(i,c)
                    xic = q_indices(i,c);
                    for kic = 1 : (Kic(i,c) - 1)
                        for kic_p = 1:Kic(i,c)
                            if kic ~= kic_p
                                rate = PH{i}{c}{1}(kic,kic_p);
                                if ni > nservers(i)
                                    dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate* nservers(i) * w(i,c)/wni;
                                    dx(xic+kic_p-1) = dx(xic+kic_p-1) + x(xic+kic-1)*rate* nservers(i) * w(i,c)/wni;
                                else
                                    dx(xic+kic-1) = dx(xic+kic-1) - x(xic+kic-1)*rate;
                                    dx(xic+kic_p-1) = dx(xic+kic_p-1) + x(xic+kic-1)*rate;
                                end
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
                                                rate = w(i,c)/wni * nservers(i) * rate;
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