classdef EmpiricalCDF < Distribution
    % Empirical Cdf for a distribution
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        data;
    end

    methods
        %Constructor
        function self = EmpiricalCDF(xdata,cdfdata)
            % SELF = EMPIRICAL(data)
            self@Distribution('EmpiricalCdf',2,[-Inf,Inf]);
            if nargin ==1
                self.data = xdata;
            else
                self.data = unique([cdfdata,xdata],'rows');
            end
        end

        function X = sample(self, n)
            % X = SAMPLE(n)
            % Get n samples from the distribution
            Y = rand(n,1);
            xset = find(self.data(:,1)<1); % remove duplicates when numerically F(x)=1
            X = interp1(self.data(xset,1),self.data(xset,2),Y(:),'spline','extrap');
        end

        function RATE = getRate(self)
            % RATE = GETRATE()
            % Get distribution rate
            RATE = 1 / getMean(self);
        end

        function MEAN = getMean(self)
            % MEAN = GETMEAN()
            % Get distribution mean
            MEAN = getMoments(self);
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            [~,~,~,SCV,~] = getMoments(self);
        end

        function VAR = getVar(self)
            % VAR = GETVAR()
            % Get distribution variance
            VAR = getSCV(self)*getMean(self)^2;
        end

        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            % Get distribution skewness
            [~,~,~,~,SKEW] = getMoments(self);
        end

        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            % Evaluate the cumulative distribution function at t
            xset = find(self.data(:,1)<1); % remove duplicates when numerically F(x)=1
            self.data(xset,1),self.data(xset,2)
            Ft = interp1(self.data(xset,1),self.data(xset,2),t,'spline','extrap');
        end

        function L = evalLST(self, s)
            % L = EVALLST(S)
            % Evaluate the Laplace transform of the distribution function at t
            cdfdata = self.data;
            row = size(cdfdata,1);
            L = 0; % the first moment
            for i = 1:1:row-1
                x = ((cdfdata(i+1,2)-cdfdata(i,2))/2+cdfdata(i,2));
                bin1 = exp(-s*x)*((cdfdata(i+1,1)-cdfdata(i,1)));
                L = L+bin1;
            end
        end

        function [m1,m2,m3,SCV,SKEW] = getMoments(self)
            cdfdata = self.data;
            row = size(cdfdata,1);
            m1 = 0; % the first moment
            m2 = 0; % the second moment
            m3 = 0; % the third moment
            for i = 1:1:row-1
                x = ((cdfdata(i+1,2)-cdfdata(i,2))/2+cdfdata(i,2));
                bin1 = x*((cdfdata(i+1,1)-cdfdata(i,1)));
                bin2 = x^2*((cdfdata(i+1,1)-cdfdata(i,1)));
                bin3 = x^3*((cdfdata(i+1,1)-cdfdata(i,1)));
                m1 = m1+bin1;
                m2 = m2+bin2;
                m3 = m3+bin3;
            end
            SCV = (m2/m1^2)-1;
            SKEW = (m3-3*m1*(m2-m1^2)-m1^3)/((m2-m1^2)^(3/2));
        end
    end
end

