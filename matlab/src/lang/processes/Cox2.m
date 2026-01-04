classdef (Sealed)  Cox2 < Markovian
    % Static class to fit two-phase coxian statistical distribution
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)
        function cx  = fit(MEAN,SCV,SKEW)
            % CX = FIT(MEAN,SCV,SKEW)
            
            cx = Cox2.fitCentral(MEAN,SCV*MEAN^2,SKEW);

            jline_mu = java.util.ArrayList();
            jline_phi = java.util.ArrayList();
            if length(line_dist.params) == 3
                jline_mu.add(cx.getParam(1).paramValue);
                jline_mu.add(cx.getParam(2).paramValue);
                jline_phi.add(cx.getParam(3).paramValue);
            else
                mu = line_dist.getParam(1).paramValue;
                phi = line_dist.getParam(2).paramValue;
                for i = 1:length(mu)
                    jline_mu.add(mu(i));
                end

                for i = 1:length(phi)
                    jline_phi.add(phi(i));
                end
            end
            cx.obj = jline.lang.processes.Coxian(jline_mu, jline_phi);
        end
        
        function cx = fitCentral(MEAN,VAR,SKEW)
            % CX = FITCENTRAL(MEAN,VAR,SKEW)
            
            % Fit the distribution from first three central moments (mean,
            % variance, skewness)
            SCV = VAR/MEAN^2;
            if nargin == 2
                cx = Cox2.fitMeanAndSCV(MEAN,SCV);
                return
            end
            e1 = MEAN;
            e2 = (1+SCV)*e1^2;
            e3 = -(2*e1^3-3*e1*e2-SKEW*(e2-e1^2)^(3/2));
            % consider the two possible solutions
            phi = (6*e1^3 - 6*e2*e1 + e3)/(- 6*e1^3 + 3*e2*e1);
            mu1 = [
                (2*(e3 - 3*e1*e2))/(- 3*e2^2 + 2*e1*e3) + (3*e1*e2 - e3 + (24*e1^3*e3 - 27*e1^2*e2^2 - 18*e1*e2*e3 + 18*e2^3 + e3^2)^(1/2))/(- 3*e2^2 + 2*e1*e3)
                (2*(e3 - 3*e1*e2))/(- 3*e2^2 + 2*e1*e3) - (e3 - 3*e1*e2 + (24*e1^3*e3 - 27*e1^2*e2^2 - 18*e1*e2*e3 + 18*e2^3 + e3^2)^(1/2))/(- 3*e2^2 + 2*e1*e3)
                ];
            mu2 = [
                -(3*e1*e2 - e3 + (24*e1^3*e3 - 27*e1^2*e2^2 - 18*e1*e2*e3 + 18*e2^3 + e3^2)^(1/2))/(- 3*e2^2 + 2*e1*e3)
                (e3 - 3*e1*e2 + (24*e1^3*e3 - 27*e1^2*e2^2 - 18*e1*e2*e3 + 18*e2^3 + e3^2)^(1/2))/(- 3*e2^2 + 2*e1*e3)];
            if phi>=0 && phi<=1 && mu1(1) >= 0 && mu2(1) >= 0
                % if the first solution is feasible
                cx = Coxian([mu1(1),mu2(1)],[phi,1]);
            elseif phi >=0 && phi <=1 && mu1(2) >= 0 && mu2(2) >= 0
                % if the second solution is feasible
                cx = Coxian([mu1(2),mu2(2)],[phi,1]);
            else
                line_warning(mfilename,'Cox2.fitCentral: Third moment could not be fitted exactly.\n');
                % fit is not feasible
                if SCV>=0.5
                    %line_warning(mfilename,'Infeasible combination of central moments, fitting only mean and squared coefficient of variation.');
                    cx = Cox2.fitMeanAndSCV(MEAN, SCV);
                else
                    %line_warning(mfilename,'Infeasible combination of central moments, fitting only mean.');
                    cx = Cox2.fitMean(MEAN);
                end
            end
        end
        
        function cx = fitMean(MEAN)
            % CX = FITMEAN(MEAN)
            
            p = 1.0-GlobalConstants.CoarseTol;
            l0 = 1/MEAN;
            l1 = 1/MEAN;
            cx = Coxian([l0,l1],[p,1]);
        end
        
        function cx = fitMeanAndSCV(MEAN, SCV)
            % CX = FITMEANANDSCV(MEAN, SCV)
            
            % Fit a 2-phase Coxian distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            if (SCV <1 && SCV >=0.5)
                l0 = 2/MEAN/(1+sqrt(1+2*(SCV-1)));
                l1 = 2/MEAN/(1-sqrt(1+2*(SCV-1)));
                p = 0.0;
            elseif (SCV == 1.0)
                l0 = 1/MEAN;
                l1 = 1/MEAN;
                p = 1.0;
            else
                l0 = 2/MEAN;
                l1 = l0/(2*SCV);
                p = 1 - l1/l0;
            end
            cx = Coxian([l0,l1],[p,1]);
        end
    end
    
end
