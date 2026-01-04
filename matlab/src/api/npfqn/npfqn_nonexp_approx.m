function [ST,gamma,nservers,rho,scva,scvs,eta] = npfqn_nonexp_approx(method,sn,ST,V,SCV,T,U,gamma,nservers)
% handler for non-exponential service and arrival processes
M = sn.nstations;
rho = zeros(M,1);
scva = ones(M,1);
scvs = ones(M,1);
eta = ones(M,1);
switch method
    case {'default','none'}
        % no-op
    case {'hvmva'}
        % no-op
    case {'interp'}
        for ist=1:M
            nnzClasses = isfinite(ST(ist,:)) & isfinite(SCV(ist,:));
            rho(ist) = sum(U(ist,nnzClasses));
            if ~isempty(nnzClasses) && any(nnzClasses)
                switch sn.sched(ist)
                    case SchedStrategy.FCFS
                        if range(ST(ist,nnzClasses))>0 || (max(SCV(ist,nnzClasses))>1 + GlobalConstants.FineTol || min(SCV(ist,nnzClasses))<1 - GlobalConstants.FineTol) % check if non-product-form
                            scva(ist) = 1; %use a M/G/k approximation
                            scvs(ist) = (SCV(ist,nnzClasses)*T(ist,nnzClasses)')/sum(T(ist,nnzClasses));
                            % multi-server asymptotic decay rate
                            gamma(ist) = (rho(ist)^nservers(ist)+rho(ist))/2;

                            if scvs(ist) > 1-1e-6 && scvs(ist) < 1+1e-6 && nservers(ist)==1
                                eta(ist) = rho(ist);
                                %continue % use M/M/1
                            else
                                % single-server (diffusion approximation, Kobayashi JACM)
                                eta(ist) = exp(-2*(1-rho(ist))/(scvs(ist)+scva(ist)*rho(ist)));
                                %[~,eta(i)]=qsys_gig1_approx_klb(sum(T(i,nnzClasses)),sum(T(i,nnzClasses))/rho(i),sqrt(scva(i)),sqrt(scvs(i)));
                            end
                            % interpolation (Sec. 4.2, LINE paper at WSC 2020)
                            % ai, bi coefficient here set to the 8th power as
                            % numerically appears to be better than 4th power
                            order = 8;
                            ai = rho(ist)^order;
                            bi = rho(ist)^order;
                            % for all classes
                            for k=find(nnzClasses)
                                if sn.rates(ist,k)>0
                                    ST(ist,k) = max(0,1-ai)*ST(ist,k) + ai*(bi*eta(ist) + max(0,1-bi)*gamma(ist))*(nservers(ist)/sum(T(ist,nnzClasses)));
                                end
                            end
                            % we are already account for multi-server effects
                            % in the scaled service times
                            nservers(ist) = 1;
                        end
                end
            end
        end
end % method
end
