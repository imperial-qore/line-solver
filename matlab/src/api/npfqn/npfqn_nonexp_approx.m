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
        for i=1:M
            nnzClasses = isfinite(ST(i,:)) & isfinite(SCV(i,:));
            rho(i) = sum(U(i,nnzClasses));
            if ~isempty(nnzClasses) && any(nnzClasses)
                switch sn.schedid(i)
                    case SchedStrategy.ID_FCFS
                        if range(ST(i,nnzClasses))>0 || (max(SCV(i,nnzClasses))>1 + GlobalConstants.FineTol || min(SCV(i,nnzClasses))<1 - GlobalConstants.FineTol) % check if non-product-form
                            scva(i) = 1; %use a M/G/k approximation
                            scvs(i) = (SCV(i,nnzClasses)*T(i,nnzClasses)')/sum(T(i,nnzClasses));
                            % multi-server asymptotic decay rate
                            gamma(i) = (rho(i)^nservers(i)+rho(i))/2;

                            if scvs(i) > 1-1e-6 && scvs(i) < 1+1e-6 && nservers(i)==1
                                eta(i) = rho(i);
                                %continue % use M/M/1
                            else
                                % single-server (diffusion approximation, Kobayashi JACM)
                                eta(i) = exp(-2*(1-rho(i))/(scvs(i)+scva(i)*rho(i)));
                                %[~,eta(i)]=qsys_gig1_approx_klb(sum(T(i,nnzClasses)),sum(T(i,nnzClasses))/rho(i),sqrt(scva(i)),sqrt(scvs(i)));
                            end
                            % interpolation (Sec. 4.2, LINE paper at WSC 2020)
                            % ai, bi coefficient here set to the 8th power as
                            % numerically appears to be better than 4th power
                            order = 8;
                            ai = rho(i)^order;
                            bi = rho(i)^order;
                            % for all classes
                            for k=find(nnzClasses)
                                if sn.rates(i,k)>0
                                    ST(i,k) = max(0,1-ai)*ST(i,k) + ai*(bi*eta(i) + max(0,1-bi)*gamma(i))*(nservers(i)/sum(T(i,nnzClasses)));
                                end
                            end
                            % we are already account for multi-server effects
                            % in the scaled service times
                            nservers(i) = 1;
                        end
                end
            end
        end
end % method
end
