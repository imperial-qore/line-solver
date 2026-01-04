classdef Distribution < Copyable
    % Distribution Abstract base class for statistical distributions
    %
    % Distribution provides the common interface and functionality for all
    % statistical distributions used in queueing models. It defines abstract
    % methods for sampling, computing moments, and obtaining distribution
    % properties that must be implemented by concrete distribution classes.
    %
    % @brief Abstract base class for all statistical distribution types
    %
    % Key characteristics:
    % - Abstract interface for all distributions
    % - Support for both discrete and continuous distributions
    % - Moment computation (mean, variance, SCV)
    % - Random sampling capabilities
    % - Parameter management and validation
    % - Copyable interface for cloning distributions
    %
    % Distribution hierarchy includes:
    % - Continuous distributions (Exp, Erlang, HyperExp, etc.)
    % - Discrete distributions (Poisson, Bernoulli, Zipf, etc.)  
    % - Phase-type distributions (PH, APH, Coxian, etc.)
    % - Markov-modulated processes (MAP, MMPP, etc.)
    %
    % Common usage patterns:
    % @code
    % dist = Exp(1.5);           % Exponential with rate 1.5
    % mean_val = dist.getMean(); % Get mean value
    % samples = dist.sample(100); % Generate 100 samples
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Hidden)
        mean; % cached mean
        immediate; % true if immediate
    end
    
    properties
        obj 
        name
        params
        support; % support interval
    end

    methods %(Abstract)
        
        function X = sample(self,n)
            % SAMPLE Generate random samples from the distribution
            %
            % @brief Generates n random samples from this distribution
            % @param n Number of samples to generate
            % @return X Vector of n random samples
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        
        function RATE = getRate(self)
            % GETRATE Get the rate parameter (inverse of mean)
            %
            % @brief Returns the rate parameter defined as 1/mean
            % @return RATE Rate parameter (1/mean)
            RATE = 1 / getMean(self);
        end
        
        function MEAN = getMean(self)
            % MEAN = GETMEAN()
            % Get distribution mean
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get distribution squared coefficient of variation (SCV = variance / mean^2)
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        
        function VAR = getVar(self)
            % VAR = GETVAR()
            % Get distribution variance
            VAR = getSCV(self)*getMean(self)^2;
        end
        
        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            % Get distribution skewness
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        
        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            % Evaluate the cumulative distribution function at t
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
        
        function L = evalLST(self, s)
            % L = EVALLAPLACETRANSFORM(S)
            % Evaluate the Laplace transform of the distribution function at t
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
            
        end
    end
    
    methods (Hidden)
        function self = Distribution(name, numParam, support)
            % SELF = DISTRIB(NAME, NUMPARAM, SUPPORT)
            % Construct a distribution from name, number of parameters, and
            % range
            self.name = name;
            self.support = support;
            self.setNumParams(numParam);
            self.immediate = false;
        end
    end
    
    methods
        function nParam = setNumParams(self, numParam)
            % NPARAM = SETNUMPARAMS(NUMPARAM)
            % Initializes the parameters
            self.params = cell(1,numParam);
            for i=1:numParam
                self.params{i}=struct('paramName','','paramValue',NaN);
            end
        end
        
        function nParam = getNumParams(self)
            % NPARAM = GETNUMPARAMS()
            % Returns the number of parameters needed to specify the distribution
            nParam = length(self.params);
        end
        
        function setParam(self, id, name, value)
            % SETPARAM(ID, NAME, VALUE, TYPECLASS)
            % Set a distribution parameter given id, name, value, Java
            % class type (for JMT translation)
            self.params{id}.paramName=name;
            self.params{id}.paramValue=value;
        end
        
        function param = getParam(self,id)
            % PARAM = GETPARAM(SELF,ID)
            % Return the parameter associated to the given id
            param = self.params{id};
        end
        
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()
            % Check if the distribution is equivalent to a Disabled
            % distribution
            %bool = cellfun(@(c) isnan(c.paramValue), self.params)
            bool = isnan(getMean(self));
        end
        
        function bool = isImmediate(self)
            % BOOL = ISIMMEDIATE()
            % Check if the distribution is equivalent to an Immediate
            % distribution
            bool = self.immediate;
        end
        
        function bool = isContinuous(self)
            % BOOL = ISCONTINUOUS()
            % Check if the distribution is discrete
            bool = isa(self,'ContinuousDistribution');
        end
        
        function bool = isDiscrete(self)
            % BOOL = ISDISCRETE()
            % Check if the distribution is discrete
            bool = isa(self,'DiscreteDistribution');
        end
        
        function delta = evalProbInterval(self,t0,t1)
            % DELTA = EVALPROBINTERVAL(SELF,T0,T1)
            % Evaluate the probability mass between t0 and t1 (t1>t0)
            if t1>=t0
                Ft1 = self.evalCDF(t1);
                Ft0 = self.evalCDF(t0);
                delta = Ft1 - Ft0;
            else
                line_error(mfilename,'CDF interval incorrectly specified (t1<t0)');
            end
        end
    end
end
