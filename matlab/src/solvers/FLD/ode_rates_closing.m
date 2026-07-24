function rates = ode_rates_closing(x, M, K, enabled, q_indices, Kic, nservers, w, sched_id, rateBase, eventIdx)
% RATES = ODE_RATES_BOTTQ(X, M, K, Q_INDICES, KIC, NSERVERS, W, STRATEGY, RATEBASE, EVENTIDX)

rates = x; % basic vector valid for INF and PS case min(ni,nservers(i))=ni
for i = 1:M
    switch sched_id(i) % source
        case SchedStrategy.INF
            % do nothing
        case SchedStrategy.EXT  %EXT
            % this is treated by a delay except that we require mass
            % conservation in the local population
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                if enabled(i,k)
                    rates(idxIni) = 1-sum(x(idxIni+1:idxEnd)); % keep total mass 1 into the source for all classes at all times, not needed for idxIni+1:idxEnd as rates is initialized equal to x
                end
            end
        case SchedStrategy.PS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );
            if ni > nservers(i) % case  min = ni handled by rates = x
                rates(idxIni:idxEnd) = x(idxIni:idxEnd)/ni * nservers(i);
            end
        case SchedStrategy.FCFS
            idxIni = q_indices(i,1);
            idxEnd = q_indices(i,K) + Kic(i,K) - 1;
            ni = sum( x(idxIni:idxEnd) );
            if ni > nservers(i) % case  min = ni handled by rates = x
                rates(idxIni:idxEnd) = x(idxIni:idxEnd)/ni * nservers(i);
            end
        case SchedStrategy.DPS %DPS
            w(i,:) = w(i,:)/sum(w(i,:));
            
            ni = mean(w(i,:));
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                if enabled(i,k)
                    ni = ni + sum( w(i,k)*x(idxIni:idxEnd) );
                end
            end
            for k=1:K
                idxIni = q_indices(i,k);
                idxEnd = q_indices(i,k) + Kic(i,k) - 1;
                if enabled(i,k)
                    rates(idxIni:idxEnd) = w(i,k)*x(idxIni:idxEnd)/ni * nservers(i); % not needed for idxIni+1:idxEnd as rates is initiliazed equal to x
                end
            end
    end
end
rates = rates(eventIdx);
rates = rateBase.*rates;
end

% THIS PART TO BE KEPT AS IT ALLOWS TO MAKE RATES STATE-DEPENDENT
% function xj = get_index(j,k)
% XJ = GET_INDEX(J,K)

%     % n is the state vector
%     % j is the queue station index
%     % k is the class index
%     % RETURNS THE INDEX of the queue-length element xi! % in the state description
%     xj = 1;
%     for z = 1 : (j-1)
%         for y = 1:K
%             xj = xj + Kic(z,y);
%         end
%     end
%     for y = 1:k-1
%         xj = xj + Kic(j,y);
%     end
%
% end
%
% function rates = ode_rates(x)
% RATES = ODE_RATES(X)

%     rates = [];
%     n = zeros(1,M); % total number of jobs in each station
%     for i = 1:M
%         for c = 1:K
%             xic = q_indices(i,c);
%             n(i) = n(i) + sum( x(xic:xic+Kic(i,c)-1 ) );
%         end
%         if S(i) == sum(N)
%             n(i) = 1;
%         end
%     end
%
%
%
%     for i = 1 : M           %transition rates for departures from any station to any other station
%         for c = 1:K         %considers only transitions from the first service phase (enough for exp servers)
%             if match(i,c)>0
%             xic = q_indices(i,c);
%             for j = 1 : M
%                 for l = 1:K
%                     if P((i-1)*K+c,(j-1)*K+l) > 0
%                     for k = 1:Kic(i,c)
%                         %pure ps + fcfs correction
%                         if x(xic+k-1) > 0 && n(i) > S(i)
%                             rates = [rates; Phi{i}{c}(k) * P((i-1)*K+c,(j-1)*K+l) * Mu{i}{c}(k) * x(xic+k-1)/n(i) * S(i);]; % f_k^{dep}
%                         elseif x(xic+k-1) > 0
%                             rates = [rates; Phi{i}{c}(k) * P((i-1)*K+c,(j-1)*K+l) * Mu{i}{c}(k) * x(xic+k-1);]; % f_k^{dep}
%                         else
%                             rates = [rates; 0;]; % f_k^{dep}
%                         end
%                     end
%                     end
%                 end
%             end
%             end
%         end
%     end
%
%     for i = 1 : M           %transition rates for "next service phase" (phases 2...)
%         for c = 1:K
%             if match(i,c)>0
%             xic = q_indices(i,c);
%             for k = 1 : (Kic(i,c) - 1)
%                 if x(xic+k-1) > 0
%                     rates = [rates; (1-Phi{i}{c}(k))*Mu{i}{c}(k)*x(xic+k-1)/n(i)];
%                 else
%                     rates = [rates; 0 ]
%                 end
%             end
%             end
%         end
%     end
%
% end
%
% function d = ode_jumps(x)
% D = ODE_JUMPS(X)

%     d = [];         %returns state changes triggered by all the events
%     for i = 1 : M   %state changes from departures in service phases 2...
%         for c = 1:K
%             if match(i,c)>0
%             xic = q_indices(i,c);
%             for j = 1 : M
%                 for l = 1:K
%                     if P((i-1)*K+c,(j-1)*K+l) > 0
%                     xjl = q_indices(j,l);
%                     for k = 1 : Kic(i,c)
%                         jump = zeros(length(x),1);
%                         jump(xic) = jump(xic) - 1; %type c in stat i completes service
%                         jump(xjl) = jump(xjl) + 1; %type c job starts in stat j
%                         d = [d jump;];
%                     end
%                     end
%                 end
%             end
%             end
%         end
%     end
%
%     for i = 1 : M   %state changes from "next service phase" transition in phases 2...
%         for c = 1:K
%             if match(i,c)>0
%             xic = q_indices(i,c);
%             for k = 1 : (Kic(i,c) - 1)
%                 jump = zeros(length(x),1);
%                 jump(xic+k-1) = jump(xic+k-1) - 1;
%                 jump(xic+k) = jump(xic+k) + 1;
%                 d = [d jump;];
%             end
%             end
%         end
%     end
% end

