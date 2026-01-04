classdef Coxian < Markovian
    % The coxian statistical distribution
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = Coxian(varargin)
            % SELF = COXIAN(VARARGIN)
            % Constructs a Coxian distribution from phase rates and
            % completion probabilities, with entry probability 1 on the
            % first phase
            self@Markovian('Coxian',2);            
            if length(varargin)==2
                mu = varargin{1};
                phi = varargin{2};
                if abs(phi(end)-1)>GlobalConstants.CoarseTol && isfinite(phi(end))
                    line_error(mfilename,sprintf('The completion probability in the last Cox state must be 1.0 but it is %0.1f', phi(end)));
                end
                setParam(self, 1, 'mu', mu);
                setParam(self, 2, 'phi', phi); % completion probability in phase 1
            elseif length(varargin)==3
                mu1 = varargin{1};
                mu2 = varargin{2};
                phi1 = varargin{3};
                setParam(self, 1, 'lambda0', mu1);
                setParam(self, 2, 'lambda1', mu2);
                setParam(self, 3, 'phi0', phi1); % completion probability in phase 1
            else
                line_error(mfilename,'Coxian accepts at most 3 parameters.');
            end
            jline_mu = java.util.ArrayList();
            jline_phi = java.util.ArrayList();
            if length(self.params) == 3
                jline_mu.add(self.getParam(1).paramValue);
                jline_mu.add(self.getParam(2).paramValue);
                jline_phi.add(self.getParam(3).paramValue);
            else
                mu = self.getParam(1).paramValue;
                phi = self.getParam(2).paramValue;
                for i = 1:length(mu)
                    jline_mu.add(mu(i));
                end

                for i = 1:length(phi)
                    jline_phi.add(phi(i));
                end
            end
            self.obj = jline.lang.processes.Coxian(jline_mu, jline_phi);
            if length(self.params) == 2
                mu = getMu(self);
                phi = getPhi(self);
                ph = {diag(-mu)+diag(mu(1:end-1).*(1-phi(1:end-1)),1),[phi.*mu,zeros(length(mu),length(mu)-1)]};
            else
                mu1 = self.getParam(1).paramValue;
                mu2 = self.getParam(2).paramValue;
                phi1 = self.getParam(3).paramValue;
                ph={[-mu1,(1-phi1)*mu1;0,-mu2],[phi1*mu1,0;mu2,0]};
            end
            self.process = ph;
            self.nPhases = length(ph{1});
        end
    end

    methods
        function phases = getNumberOfPhases(self)
            % PHASES = GETNUMBEROFPHASES()
            % Return number of phases in the distribution
            if length(self.params) == 2
                phases  = length(self.getParam(1).paramValue);
            else
                phases  = 2;
            end
        end

        function ex = getMean(self)
            % EX = GETMEAN()
            % Get distribution mean
            if length(self.params) == 2
                mu = getMu(self);
                phi = getPhi(self);
                ex = map_mean({diag(-mu)+diag(mu(1:end-1).*(1-phi(1:end-1)),1),[phi.*mu,zeros(length(mu),length(mu)-1)]});
            else
                % Get distribution mean
                mu1 = self.getParam(1).paramValue;
                mu2 = self.getParam(2).paramValue;
                phi1 = self.getParam(3).paramValue;
                ex = 1/mu1 + (1-phi1)/mu2;
            end
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get the squared coefficient of variation of the distribution (SCV = variance / mean^2)
            if length(self.params) == 2
                mu = getMu(self);
                phi = getPhi(self);
                SCV = map_scv({diag(-mu)+diag(mu(1:end-1).*(1-phi(1:end-1)),1),[phi.*mu,zeros(length(mu),length(mu)-1)]});
            else
                mu1 = self.getParam(1).paramValue;
                mu2 = self.getParam(2).paramValue;
                phi1 = self.getParam(3).paramValue;
                mean = 1/mu1 + (1-phi1)/mu2;
                var = ((2*mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (2*mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1))/(mu1*mu1*((mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1))) - (1/mu1 - (phi1 - 1)/mu2)*(1/mu1 - (phi1 - 1)/mu2) - (((phi1 - 1)/(mu2*mu2) + (phi1 - 1)/(mu1*mu2))*((2*mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (2*mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1)))/((mu2*(mu1 - mu1*phi1))/(mu1 + mu2 - mu1*phi1) + (mu1*mu2*phi1)/(mu1 + mu2 - mu1*phi1));
                SCV = var / mean^2;
            end
        end

        function mu = getMu(self)
            % MU = GETMU()
            % Get vector of rates
            if length(self.params) == 2
                mu = self.getParam(1).paramValue(:);
            else
                mu1 = self.getParam(1).paramValue;
                mu2 = self.getParam(2).paramValue;
                mu = [mu1;mu2];
            end
        end

        function phi = getPhi(self)
            % PHI = GETPHI()
            % Get vector of completion probabilities
            if length(self.params) == 2
                phi = self.getParam(2).paramValue(:);
            else
                phi1 = self.getParam(3).paramValue;
                phi = [phi1;1.0];
            end
        end

    end

    methods(Static)
        function cx = fitCentral(MEAN, VAR, SKEW)
            % CX = FITCENTRAL(MEAN, VAR, SKEW)
            cx = Cox2.fitCentral(MEAN, VAR, SKEW);
            SCV = VAR/MEAN^2;
            if abs(1-map_scv(cx.getProcess)/SCV) > 0.01
                cx = Coxian.fitMeanAndSCV(MEAN, SCV);
            end
            cx.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function [cx,mu,phi] = fitMeanAndSCV(MEAN, SCV)
            % [CX,MU,PHI] = FITMEANANDSCV(MEAN, SCV)
            % Fit a Coxian distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            if SCV >= 1-GlobalConstants.CoarseTol && SCV <= 1+GlobalConstants.CoarseTol
                n = 1;
                mu = 1/MEAN;
                phi = 1;
            elseif SCV > 0.5+GlobalConstants.CoarseTol && SCV<1-GlobalConstants.CoarseTol
                phi = 0.0;
                n = 2;
                mu = zeros(n,1);
                mu(1) = 2/MEAN/(1+sqrt(1+2*(SCV-1)));
                mu(2) = 2/MEAN/(1-sqrt(1+2*(SCV-1)));
            elseif SCV <= 0.5+GlobalConstants.CoarseTol
                n = ceil(1/SCV);
                lambda = n/MEAN;
                mu = lambda*ones(n,1);
                phi = zeros(n,1);
            else % SCV > 1+GlobalConstants.CoarseTol
                n = 2;
                %transform hyperexp into coxian
                mu = zeros(n,1);
                mu(1) = 2/MEAN;
                mu(2) = mu(1)/( 2*SCV );
                phi = zeros(n,1);
                phi(1) = 1-mu(2)/mu(1);
                phi(2) = 1;
            end
            phi(n) = 1;
            cx = Coxian(mu,phi);
            cx.immediate = MEAN < GlobalConstants.CoarseTol;
        end

    end

end

