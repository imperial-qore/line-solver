classdef FunctionTask < Task 
    % A software server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        
        SetupTime; 
        SetupTimeMean;
        SetupTimeSCV;
        DelayOffTime;
        DelayOffTimeMean;
        DelayOffTimeSCV;

    end

    methods
        %public methods, including constructor
        
        %constructor
        function self = FunctionTask(model, name, multiplicity, scheduling)
            %self = FunctionTask(model, name, setupTime, delayoffTime, multiplicity, scheduling)

            if ~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end
                       
            if nargin < 5
                multiplicity = 1;
            end
            
            if nargin < 6
                scheduling = SchedStrategy.FCFS;
            end    
            self@Task(model, name, multiplicity, scheduling);           
            self.setSetupTime(Immediate());
            self.setDelayOffTime(Immediate());
        end

        function self = setSetupTime(self, SetupTime)
            % self = SETTHINKTIME(self, THINKTIME)
            
            if isnumeric(SetupTime)
                if SetupTime <= GlobalConstants.FineTol
                    self.SetupTime = Immediate.getInstance();
                    self.SetupTimeMean = GlobalConstants.FineTol;
                    self.SetupTimeSCV = GlobalConstants.FineTol;
                else
                    self.SetupTime = Exp(1/SetupTime);
                    self.SetupTimeMean = SetupTime;
                    self.SetupTimeSCV = 1.0;
                end
            elseif isa(SetupTime,'Distribution')
                self.SetupTime = SetupTime;
                self.SetupTimeMean = SetupTime.getMean();
                self.SetupTimeSCV = SetupTime.getSCV();
            end
        end

        function self = setDelayOffTime(self, DelayOffTime)
            % self = SETTHINKTIME(self, THINKTIME)
            
            if isnumeric(DelayOffTime)
                if DelayOffTime <= GlobalConstants.FineTol
                    self.DelayOffTime = Immediate.getInstance();
                    self.DelayOffTimeMean = GlobalConstants.FineTol;
                    self.DelayOffTimeSCV = GlobalConstants.FineTol;
                else
                    self.DelayOffTime = Exp(1/DelayOffTime);
                    self.DelayOffTimeMean = DelayOffTime;
                    self.DelayOffTimeSCV = 1.0;
                end
            elseif isa(DelayOffTime,'Distribution')
                self.DelayOffTime = DelayOffTime;
                self.DelayOffTimeMean = DelayOffTime.getMean();
                self.DelayOffTimeSCV = DelayOffTime.getSCV();
            end
        end
    end
end