%{
%{
 % @file pfqn_stdf_heur.m
 % @brief Heuristic sojourn time distribution for multiserver FCFS nodes (McKenna 1987 variant).
%}
%}

%{
%{
 % @brief Heuristic sojourn time distribution for multiserver FCFS nodes (McKenna 1987 variant).
 % @fn pfqn_stdf_heur(L, N, Z, S, fcfsNodes, rates, tset)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param S Number of servers per station.
 % @param fcfsNodes Indices of FCFS nodes to analyze.
 % @param rates Service rates (MxR matrix).
 % @param tset Time points for CDF evaluation.
 % @return RD Cell array of sojourn time distributions {station,class}.
%}
%}
function  RD = pfqn_stdf_heur(L,N,Z,S,fcfsNodes,rates,tset)
% Heuristic sojourn time distribution analysis at multiserver FCFS nodes
% based on a variant of the method in J. McKenna 1987 JACM
stabilityWarnIssued = false;
[M,R]=size(L);
T=length(tset);
mu = zeros(M,sum(N));
for k=1:M
    for n=1:sum(N)
        mu(k,n) = min(S(k),n);
    end
end

% when t==0 the hkc function is not well defined
tset(tset == 0) = GlobalConstants.FineTol;
hkc = cell(1,R);
for k=fcfsNodes(:)'
    for r=1:R
        if L(k,r) > GlobalConstants.FineTol
            Nr = oner(N,r);
            if M==1
                options = SolverNC.defaultOptions;
                [~,lGr] = pfqn_comomrm_ld(L,Nr,Z,mu,options);
                m = 2;  [~,lGra] = pfqn_comomrm_ld(L,Nr,Z,pfqn_mu_ms(sum(N),m,S(k)),options);
                Q1 = zeros(1,R);
                for s=1:R
                    Q1(k,s) = L(k,s)* exp(lGra-lGr);
                end
                isNumStable = true;
            else
                [~,Q1,~,~,lGr,isNumStable]  = pfqn_mvald(L,Nr,Z,mu);
            end
            
            hMAPr = cell(1,1+sum(N));
            for n=0:sum(N)
                if n < S(k)
                    % the heuristic here is to use the average rates                    
                    hMAPr{1+n} = map_exponential(1/rates(k,r));
                else
                    % the heuristic here is to use the rates per class
                    C = {map_exponential(1/rates(k,r))};
                    for s=1:R
                        C = {C{:}, map_exponential(Q1(k,s)*(n-S(k)+1)/sum(Q1(k,:))/rates(k,s))};                        
                    end
                    hMAPr{1+n} = map_sumind(C);
                end
                hkc{r}(1:T,1+n) = map_cdf(hMAPr{1+n}, tset)';
            end
            
            lGr = lGr(end);
            if ~isNumStable && ~stabilityWarnIssued
                stabilityWarnIssued = true;
                line_warning(mfilename,'The computation of the sojourn time distribution is numerically unstable.\n');
            end
            RD{k,r} = zeros(length(tset),2);
            RD{k,r}(:,2) = tset(:);
            
            %% this is the original code in the paper
            %             Gkrt = zeros(T,1);
            %             nvec = pprod(Nr);
            %             while nvec >= 0
            %                 [~,~,~,~,lGk,isNumStable]  = pfqn_mvald(L(setdiff(1:M,k),:),Nr-nvec,Z,mu(setdiff(1:M,k),:));
            %                 lGk = lGk(end);
            %                 if ~isNumStable && ~stabilityWarnIssued
            %                     stabilityWarnIssued = true;
            %                     line_warning(mfilename,'The computation of the sojourn time distribution is numerically unstable');
            %                 end
            %                 for t=1:T
            %                     if sum(nvec)==0
            %                         Gkrt(t) = Gkrt(t) + hkc(t,1+0) * exp(lGk);
            %                     else
            %                         lFk = nvec*log(L(k,:))' +factln(sum(nvec)) - sum(factln(nvec)) - sum(log(mu(k,1:sum(nvec))));
            %                         Gkrt(t) = Gkrt(t) + hkc(t,1+sum(nvec)) * exp(lFk + lGk);
            %                     end
            %                 end % nvec
            %                 nvec = pprod(nvec, Nr);
            %             end
            %             lGkrt = log(Gkrt);
            %             RD{k,r}(1:T,1) = exp(lGkrt - lGr);
            
            %% this is faster as it uses the recursive form for LD models
            Hkrt = [];
            [~,~,~,~,lGk]  = pfqn_mvald(L(setdiff(1:M,k),:),Nr,Z,mu(setdiff(1:M,k),:));
            lGk = lGk(end);
            for t=1:T
                %[t,T]
                gammat = mu;
                for m=1:sum(Nr)
                    gammat(k,m) = mu(k,m) * hkc{r}(t,1+m-1) / hkc{r}(t,1+m);
                end
                gammak = mushift(gammat,k);
                Hkrt(t) = hkc{r}(t,1+0) * exp(lGk); % nvec = 0
                for s=1:R % nvec >= 1s
                    if Nr(s)>0
                        lYks_t  = pfqn_rd(L,oner(Nr,s),Z,gammak);
                        %[~,~,~,~,lYks_t]  = pfqn_mvald(L,oner(Nr,s),Z,gammak);
                        lYks_t = lYks_t(end);
                        Hkrt(t) = Hkrt(t) + (L(k,s) * hkc{r}(t,1+0) / gammat(k,1)) * exp(lYks_t); % nvec = 0
                    end
                end
            end
            Hkrt(isnan(Hkrt)) = GlobalConstants.FineTol;
            lHkrt = log(Hkrt);
            
            RD{k,r}(1:T,1) = exp(lHkrt - lGr);
        end
    end
end
end

function Fmx=Fm(m,x)
if m==1
    Fmx = 1-exp(-x);
else
    A = 0;
    for j=0:(m-1)
        A = A + x^j/factorial(j);
    end
    Fmx = 1-exp(-x)*A;
end
end

function mushifted = mushift(mu,i)
% shifts the service rate vector
[M,N]=size(mu);
for m=1:M
    if m==i
        mushifted(m,1:(N-1))=mu(m,2:N);
    else
        mushifted(m,1:(N-1))=mu(m,1:(N-1));
    end
end
end