classdef Replayer < TimeSeries
    % Empirical time series from a trace
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        data;
    end
    
    methods
        %Constructor
        function self = Replayer(data)            
            self@TimeSeries('Replayer',1);
            if ischar(data) % interpret as string
                % SELF = REPLAYER(FILENAME)
                fileName = data;
                if exist(fileName,'file') == 2
                    dirStruct = dir(fileName);
                    fileName = [dirStruct.folder,filesep,dirStruct.name];
                else
                    line_error(mfilename,'The file cannot be located, use the full file path.');
                end                
                setParam(self, 1, 'fileName', fileName);
                self.data = [];
            else
                self.data = data;
            end
        end
        
        function load(self)
            % LOAD()
            
            fileName = self.getParam(1).paramValue;
            self.data = load(fileName);
            self.data = self.data(:);
        end
        
        function unload(self)
            % UNLOAD()
            
            self.data = [];
        end
        
        function rate = getRate(self)
            rate = 1 / self.getMean;
        end
        
        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Get distribution mean
            if isempty(self.data)
                load(self);
            end
            ex = mean(self.data);
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            if isempty(self.data)
                load(self);
            end
            SCV = var(self.data)/mean(self.data)^2;
        end
        
        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            
            % Get distribution skewness
            if isempty(self.data)
                load(self);
            end
            SKEW = skewness(self.data);
        end
        
        function distr = fitExp(self)
            % DISTR = FITEXP()
            
            distr = Exp.fitMean(self.getMean);
        end
        
        function distr = fitAPH(self)
            % DISTR = FITAPH()
            distr = APH.fit(self.getMean, self.getSCV, self.getSkewness);
        end
        
        function distr = fitCoxian(self)
            % DISTR = FITCOXIAN()
            distr = Cox2.fit(self.getMean, self.getSCV, self.getSkewness);
            
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            if isempty(self.data)
                load(self);
            end
            L = mean(exp(-s*self.data));
        end
    end
end

