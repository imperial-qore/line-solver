classdef TimeSeries < PointProcess
    % An abstract class for a point process realization (a time series)
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    methods (Hidden)
        function self = TimeSeries(className, numPar)
            % SELF = TIMESERIES(CLASSNAME, NUMPAR)
            
            self@PointProcess(className,numPar);
        end
    end
    
    methods
        function transform(self, filterType, filterParam)
            % TRANSFORM(FILTERTYPE, FILTERPARAM)
            
            if isempty(self.data)
                load(self);
            end
            I = length(self.data);
            switch filterType
                case TimeSeriesFilter.Shuffle
                    self.data=self.data(randperm(I));
                case TimeSeriesFilter.MovingAvg
                    wndsz = filterParam;
                    inData = [self.data; zeros(wndsz-1,1)];
                    for i=1:I
                        self.data(i) = mean(inData(i:(i+wndsz-1)));
                    end
            end
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Returns the mean of the inter-arrival times
            ex = mean(self.data);
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get squared coefficient of variation of the interarrival times (SCV = variance / mean^2)
            SCV = var(self.data)/mean(self.data)^2;
        end
        
        function VAR = getVariance(self)
            % VAR = GETVARIANCE()
            
            % Get inter-arrival time distribution variance
            VAR = var(self.data);
        end
        
        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            
            SKEW = skewness(self.data);
            % Get skewness of the interarrival times
        end
        
        function summary(self)
            % SUMMARY()
            
            if isempty(self.data)
                load(self);
            end
            MEAN = self.getMean;
            MED = median(self.data);
            SCV = self.getSCV;
            SKEW = skewness(self.data);
            QUART = [prctile(self.data,25),prctile(self.data,50),prctile(self.data,75)];
            TAILS1PERC = [prctile(self.data,95),prctile(self.data,(1-GlobalConstants.FineTol)*100)];
            MINMAX = [min(self.data),max(self.data)];
            MAD = mad(self.data,1); %median based mad
            line_printf('\nReplayer: length=%d NaNs=%d\nMoments: mean=%f scv=%f cv=%f skew=%f\nPercentiles: p25=%f p50=%f p75=%f p95=%f\nOrder: min=%f max=%f median=%f mad=%f',length(self.data),sum(isnan(self.data)),MEAN,SCV,sqrt(SCV),SKEW,QUART(1),QUART(2),QUART(3),TAILS1PERC(1),MINMAX(1),MINMAX(2),MED,MAD);
        end
    end
end
