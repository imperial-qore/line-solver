classdef Replayer < Distribution
    % Empirical time series from a trace
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        data;
        cursample; 
    end
    
    methods
        %Constructor
        function self = Replayer(data)
            self@Distribution('Replayer',1,[0,Inf]);
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
                % Create Java object after file validation with validated path
                self.obj = jline.lang.processes.Replayer(fileName);
            else
                self.data = data;
                self.cursample = 0;
                % Create Java object directly with array data
                self.obj = jline.lang.processes.Replayer(data);
            end
        end
        
        function load(self)
            % LOAD()
            
            fileName = self.getParam(1).paramValue;
            self.data = load(fileName);
            self.data = self.data(:);
            self.cursample = 0;
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
                self.load();
            end
            ex = mean(self.data);
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            if isempty(self.data)
                self.load();
            end
            SCV = var(self.data)/mean(self.data)^2;
        end
        
        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            
            % Get distribution skewness
            if isempty(self.data)
                load(self);
            end
            SKEW = skewness(self.data,0); % ,0 ensures Apache Commons equivalence in Java
        end
        
        function distr = fitExp(self)
            % DISTR = FITEXP()
            
            distr = Exp.fitMean(self.getMean);
        end
        
        function distr = fitAPH(self)
            % DISTR = FITAPH()
            distr = APH.fit(self.getMean, self.getSCV, self.getSkewness);
                distr.obj = self.obj.fitAPH();
        end
        
        function distr = fitCoxian(self)
            % DISTR = FITCOXIAN()
            distr = Cox2.fit(self.getMean, self.getSCV, self.getSkewness);
            
        end

        function X = sample(self)
            % X = SAMPLE()
            self.cursample = self.cursample + 1;
            X = self.data(self.cursample);
        end
        
        function L = evalLST(self, s)
            % L = EVALST(S)
            % Evaluate the Laplace-Stieltjes transform of the distribution function at t
            if isempty(self.data)
                self.load();
            end
            L = mean(exp(-s*self.data));
        end
    end
end

