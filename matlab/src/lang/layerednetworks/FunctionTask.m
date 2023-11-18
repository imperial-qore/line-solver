classdef FunctionTask < Task 
    % A software server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        
        SetupTime; 
        SetupTimeMean;
        SetupTimeSCV;
        DelayedOffTime;
        DelayedOffTimeMean;
        DelayedOffTimeSCV;

    end

    methods
        %public methods, including constructor
        
        %constructor
        function self = FunctionTask(model, name, multiplicity, scheduling)
            %self = FunctionTask(model, name, setupTime, delayedoffTime, multiplicity, scheduling)

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
            self.setDelayedOffTime(Immediate());
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

        function self = setDelayedOffTime(self, DelayedOffTime)
            % self = SETTHINKTIME(self, THINKTIME)
            
            if isnumeric(DelayedOffTime)
                if DelayedOffTime <= GlobalConstants.FineTol
                    self.DelayedOffTime = Immediate.getInstance();
                    self.DelayedOffTimeMean = GlobalConstants.FineTol;
                    self.DelayedOffTimeSCV = GlobalConstants.FineTol;
                else
                    self.DelayedOffTime = Exp(1/DelayedOffTime);
                    self.DelayedOffTimeMean = DelayedOffTime;
                    self.DelayedOffTimeSCV = 1.0;
                end
            elseif isa(DelayedOffTime,'Distribution')
                self.DelayedOffTime = DelayedOffTime;
                self.DelayedOffTimeMean = DelayedOffTime.getMean();
                self.DelayedOffTimeSCV = DelayedOffTime.getSCV();
            end
        end
    end
end