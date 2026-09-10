classdef Logger < Node
    % A node where jobs are logged upon passage.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        fileName;
        filePath;
        schedPolicy;
        schedStrategy;
        cap;
    end
    
    properties (Access=private)
        wantStartTime;
        wantLoggerName;
        wantTimestamp;
        wantJobID;
        wantJobClass;
        wantTimeSameClass;
        wantTimeAnyClass;
    end
    
    methods
        %Constructor
        function self = Logger(model, name, logFileName)
            % SELF = LOGGER(MODEL, NAME, LOGFILENAME)
            
            self@Node(name);
            if model.isMatlabNative()
                [~,fileName,fileExt] = fileparts(logFileName);
                self.fileName = sprintf('%s%s',fileName,fileExt);
                if isempty(model.getLogPath)
                    line_error(mfilename,'To instantiate a Logger, first use setLogPath method on the Network object to define the global path to save logs.');
                else
                    self.filePath = model.getLogPath;
                end
                classes = model.getClasses();
                self.input = Buffer(classes);
                self.output = Dispatcher(classes);
                self.cap = Inf;
                self.schedPolicy = SchedStrategyType.NP;
                self.schedStrategy = SchedStrategy.FCFS;
                self.server = LogTunnel();
                self.setStartTime(false);
                self.setLoggerName(false);
                self.setTimestamp(true);
                self.setJobID(true);
                self.setJobClass(true);
                self.setTimeSameClass(false);
                self.setTimeAnyClass(false);
                self.setModel(model);
                self.model.addNode(self);
            elseif model.isJavaNative()
                self.setModel(model);
                self.obj = jline.lang.nodes.Logger(model.obj, name, logFileName);
                self.index = model.obj.getNodeIndex(self.obj);
            end
        end
        
        function ret = getStartTime(self)
            % RET = GETSTARTTIME()
            
            ret = self.wantStartTime;
        end
        function ret = getLoggerName(self)
            % RET = GETLOGGERNAME()
            
            ret = self.wantLoggerName;
        end
        function ret = getTimestamp(self)
            % RET = GETTIMESTAMP()
            
            ret = self.wantTimestamp;
        end
        function ret = getJobID(self)
            % RET = GETJOBID()
            
            ret = self.wantJobID;
        end
        function ret = getJobClass(self)
            % RET = GETJOBCLASS()
            
            ret = self.wantJobClass;
        end
        function ret = getTimeSameClass(self)
            % RET = GETTIMESAMECLASS()
            
            ret = self.wantTimeSameClass;
        end
        function ret = getTimeAnyClass(self)
            % RET = GETTIMEANYCLASS()
            
            ret = self.wantTimeAnyClass;
        end
        
        function setStartTime(self, bool)
            % SETSTARTTIME(BOOL)
            
            if bool
                self.wantStartTime = 'true';
            else
                self.wantStartTime = 'false';
            end
        end
        
        function setTimestamp(self, bool)
            % SETTIMESTAMP(BOOL)
            
            if bool
                self.wantTimestamp = 'true';
            else
                self.wantTimestamp = 'false';
            end
        end
        
        function setLoggerName(self, bool)
            % SETLOGGERNAME(BOOL)
            
            if bool
                self.wantLoggerName = 'true';
            else
                self.wantLoggerName = 'false';
            end
        end
        
        function setTimeSameClass(self, bool)
            % SETTIMESAMECLASS(BOOL)
            
            if bool
                self.wantTimeSameClass = 'true';
            else
                self.wantTimeSameClass = 'false';
            end
        end
        
        function setTimeAnyClass(self, bool)
            % SETTIMEANYCLASS(BOOL)
            
            if bool
                self.wantTimeAnyClass = 'true';
            else
                self.wantTimeAnyClass = 'false';
            end
        end
        
        function setJobID(self, bool)
            % SETJOBID(BOOL)
            
            if bool
                self.wantJobID = 'true';
            else
                self.wantJobID = 'false';
            end
        end
        
        function setJobClass(self, bool)
            % SETJOBCLASS(BOOL)
            
            if bool
                self.wantJobClass = 'true';
            else
                self.wantJobClass = 'false';
            end
        end
        
        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)
            
            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end
        
    end
end
