classdef SolverLQNS < LayeredNetworkSolver
    % A solver that interfaces the LQNS to LINE.
    %
    % Copyright (c) 2012-2022, Imperial College London
    % All rights reserved.

    % Taken from: https://gitlab.doc.ic.ac.uk/zn19/line-dev/-/blob/master/src/solvers/@SolverLQNS/SolverLQNS.m?ref_type=heads
    
    methods
        function self = SolverLQNS(model, varargin)
            % SELF = SOLVERLQNS(MODEL, VARARGIN)
            
            self@LayeredNetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            if ~SolverLQNS.isAvailable()
                line_error(mfilename,'SolverLQNS requires the lqns and lqsim commands to be available on the system path. Please visit: http://www.sce.carleton.ca/rads/lqns/');
            end
        end
%---------------------------------------------------------------------------------------------------        
        function runtime = runAnalyzer(self, options)
            % RUNTIME = RUN()
            % Run the solver
            
            tic;
            options = self.getOptions;
            filename = [lineTempName,'.lqnx'];
            self.model.writeXML(filename);
            if options.verbose
                ignoreWarn = '';
            else
                ignoreWarn = '-w ';
            end
            
            if isunix
                switch options.method
                    case {'default','lqns'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvn'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Playering=srvn -Pstop-on-message-loss=false -x ',filename]);
                    case {'exact'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvnexact'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'sim','lqsim'}
                        system(['lqsim ',ignoreWarn,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename]);
                    case {'lqnsdefault'}
                        system(['lqns ',ignoreWarn,' -x ',filename]);
                    otherwise
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                end
            else
                switch options.method
                    case {'default','lqns'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvn'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Playering=srvn -Pstop-on-message-loss=false -x ',filename]);
                    case {'exact'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvnexact'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'sim','lqsim'}
                        system(['lqsim ',ignoreWarn,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename]);
                    case {'lqnsdefault'}
                        system(['lqns ',ignoreWarn,' -x ',filename]);
                    otherwise
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                end
            end
            
            self.parseXMLResults(filename);
            
            if ~options.keep
                [filepath,name] = fileparts(filename);
                delete([filepath,filesep,name,'.*'])
            end
            runtime = toc;
        end
%----------------------------------------------------------------------------------------------------        
        function [runtime,exceedProb] = runAnalyzer_thresholdTime(self, thresholdTime)
            % RUNTIME = RUN()
            % Run the solver
            
            tic;
            options = self.getOptions;
            filename = [lineTempName,'.lqnx'];
%             self.model.writeXML(filename);
            self.model.writeXMLthresholdTime(filename, 0, thresholdTime);
            if options.verbose
                ignoreWarn = '';
            else
                ignoreWarn = '-w ';
            end
            
            if isunix
                switch options.method
                    case {'default','lqns'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvn'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Playering=srvn -Pstop-on-message-loss=false -x ',filename]);
                    case {'exact'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvnexact'}
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'sim','lqsim'}
                        system(['lqsim ',ignoreWarn,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename]);
                    case {'lqnsdefault'}
                        system(['lqns ',ignoreWarn,' -x ',filename]);
                    otherwise
                        system(['lqns ',ignoreWarn,' --iteration-limit ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                end
            else
                switch options.method
                    case {'default','lqns'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvn'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Playering=srvn -Pstop-on-message-loss=false -x ',filename]);
                    case {'exact'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'srvnexact'}
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Playering=srvn -Pmva=exact -Pstop-on-message-loss=false -x ',filename]);
                    case {'sim','lqsim'}
                        system(['lqsim ',ignoreWarn,' -A ',num2str(options.samples),',3  -Pstop-on-message-loss=false -x ',filename]);
                    case {'lqnsdefault'}
                        system(['lqns ',ignoreWarn,' -x ',filename]);
                    otherwise
                        system(['lqns ',ignoreWarn,' -i ',num2str(options.iter_max),' -Pstop-on-message-loss=false -x ',filename]);
                end
            end
            
%             self.parseXMLResults(filename);
            [~,~,exceedProb] = self.parseXMLResultsExceedProb(filename);
            
            if ~options.keep
                [filepath,name] = fileparts(filename);
                delete([filepath,filesep,name,'.*'])
            end
            runtime = toc;
        end
%----------------------------------------------------------------------------------------------------        
        function [QN,UN,RN,SCVN,TN,AN,WN,runningtime] = getAvg(self,~,~,~,~,useLQNSnaming)
            % [QN,UN,RN,TN] = GETAVG(SELF,~,~,~,~,USELQNSnaming)
            %
            % SN: average service time
            
            if nargin < 5
                useLQNSnaming = false;
            end
            
            runningtime = runAnalyzer(self);
            QN = self.result.Avg.QLen;
            UN = self.result.Avg.Util;
            RN = self.result.Avg.RespT;
            SCVN = self.result.Avg.SCV;
            TN = self.result.Avg.Tput;
            PN = self.result.Avg.ProcUtil;
            SN = self.result.Avg.SvcT;
            AN = TN*0;
            WN = RN*0;
            
            if ~useLQNSnaming
                QN = UN;
                UN = PN;
                RN = SN;
            end
            
        end
 %------------------------------------------------------------------------------------------------
        function [result, iterations] = parseXMLResults(self, filename)
            % [RESULT, ITERATIONS] = PARSEXMLRESULTS(FILENAME)
            
            import javax.xml.parsers.*;
            import org.w3c.dom.*;
            import java.io.*;
            
            lqn = self.getStruct;
            numOfNodes = lqn.nidx;
            numOfCalls = lqn.ncalls;
            Avg.Nodes.Utilization = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase1Utilization = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase2Utilization = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase1ServiceTime = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase1ServiceSCV = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase2ServiceTime = NaN*ones(numOfNodes,1);
            Avg.Nodes.Throughput = NaN*ones(numOfNodes,1);
            Avg.Nodes.ProcWaiting = NaN*ones(numOfNodes,1);
            Avg.Nodes.ProcUtilization = NaN*ones(numOfNodes,1);
            Avg.Edges.Waiting = NaN*ones(numOfCalls,1);
            
            % init Java XML parser and load file
            dbFactory = DocumentBuilderFactory.newInstance();
            dBuilder = dbFactory.newDocumentBuilder();
            
            [fpath,fname,~] = fileparts(filename);
            resultFilename = [fpath,filesep,fname,'.lqxo'];
            if self.options.verbose
                line_printf('\nParsing LQNS result file: %s',resultFilename);
            end
            if self.options.keep
                line_printf('\nLQNS result file available at: %s',resultFilename);
            end
            
            doc = dBuilder.parse(resultFilename);
            doc.getDocumentElement().normalize();
            
            %solver-params
            solverParams = doc.getElementsByTagName('solver-params');
            for i = 0:solverParams.getLength()-1
                solverParam = solverParams.item(i);
                result = solverParam.getElementsByTagName('result-general');
                iterations = str2double(result.item(0).getAttribute('iterations'));
            end
            
            procList = doc.getElementsByTagName('processor');
            for i = 0:procList.getLength()-1
                %Element - Host
                procElement = procList.item(i);
                procName = char(procElement.getAttribute('name'));
                procPos = findstring(lqn.names,procName);
                procResult = procElement.getElementsByTagName('result-processor');
                uRes = str2double(procResult.item(0).getAttribute('utilization'));
                Avg.Nodes.ProcUtilization(procPos) = uRes;
                
                taskList = procElement.getElementsByTagName('task');
                for j = 0:taskList.getLength()-1
                    %Element - Task
                    taskElement = taskList.item(j);
                    taskName = char(taskElement.getAttribute('name'));
                    taskPos = findstring(lqn.names,taskName);
                    taskResult = taskElement.getElementsByTagName('result-task');
                    uRes = str2double(taskResult.item(0).getAttribute('utilization'));
                    p1uRes = str2double(taskResult.item(0).getAttribute('phase1-utilization'));
                    p2uRes = str2double(taskResult.item(0).getAttribute('phase2-utilization'));
                    tRes = str2double(taskResult.item(0).getAttribute('throughput'));
                    puRes = str2double(taskResult.item(0).getAttribute('proc-utilization'));
                    Avg.Nodes.Utilization(taskPos) = uRes;
                    Avg.Nodes.Phase1Utilization(taskPos) = p1uRes;
                    Avg.Nodes.Phase2Utilization(taskPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                    Avg.Nodes.Throughput(taskPos) = tRes;
                    Avg.Nodes.ProcUtilization(taskPos) = puRes;
                    
                    entryList = taskElement.getElementsByTagName('entry');
                    for k = 0:entryList.getLength()-1
                        %Element - Entry
                        entryElement = entryList.item(k);
                        entryName = char(entryElement.getAttribute('name'));
                        entryPos = findstring(lqn.names,entryName);
                        entryResult = entryElement.getElementsByTagName('result-entry');
                        uRes = str2double(entryResult.item(0).getAttribute('utilization'));
                        p1uRes = str2double(entryResult.item(0).getAttribute('phase1-utilization'));
                        p1SCV = str2double(entryResult.item(0).getAttribute('squared-coeff-variation'));
                        p2uRes = str2double(entryResult.item(0).getAttribute('phase2-utilization'));
                        p1stRes = str2double(entryResult.item(0).getAttribute('phase1-service-time'));
                        p2stRes = str2double(entryResult.item(0).getAttribute('phase2-service-time'));
                        tRes = str2double(entryResult.item(0).getAttribute('throughput'));
                        puRes = str2double(entryResult.item(0).getAttribute('proc-utilization'));
                        Avg.Nodes.Utilization(entryPos) = uRes;
                        Avg.Nodes.Phase1Utilization(entryPos) = p1uRes;
                        Avg.Nodes.Phase2Utilization(entryPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                        Avg.Nodes.Phase1ServiceTime(entryPos) = p1stRes;
                        Avg.Nodes.Phase1ServiceSCV(entryPos) = p1SCV;
                        Avg.Nodes.Phase2ServiceTime(entryPos) = ifthenelse(isempty(p2stRes),NaN,p2stRes);
                        Avg.Nodes.Throughput(entryPos) = tRes;
                        Avg.Nodes.ProcUtilization(entryPos) = puRes;
                    end
                    
                    %task-activities
                    taskActsList = taskElement.getElementsByTagName('task-activities');
                    if taskActsList.getLength > 0
                        taskActsElement = taskActsList.item(0);
                        actList = taskActsElement.getElementsByTagName('activity');
                        for l = 0:actList.getLength()-1
                            %Element - Activity
                            actElement = actList.item(l);
                            if strcmp(char(actElement.getParentNode().getNodeName()),'task-activities')
                                actName = char(actElement.getAttribute('name'));
                                actPos = findstring(lqn.names,actName);
                                actResult = actElement.getElementsByTagName('result-activity');
                                uRes = str2double(actResult.item(0).getAttribute('utilization'));
                                stRes = str2double(actResult.item(0).getAttribute('service-time'));
                                tRes = str2double(actResult.item(0).getAttribute('throughput'));
                                pwRes = str2double(actResult.item(0).getAttribute('proc-waiting'));
                                puRes = str2double(actResult.item(0).getAttribute('proc-utilization'));
                                Avg.Nodes.Utilization(actPos) = uRes;
                                Avg.Nodes.Phase1ServiceTime(actPos) = stRes;
                                Avg.Nodes.Throughput(actPos) = tRes;
                                Avg.Nodes.ProcWaiting(actPos) = pwRes;
                                Avg.Nodes.ProcUtilization(actPos) = puRes;
                                
                                actID = lqn.names{actPos};
                                %synch-call
                                synchCalls = actElement.getElementsByTagName('synch-call');
                                for m = 0:synchCalls.getLength()-1
                                    callElement = synchCalls.item(m);
                                    destName = char(callElement.getAttribute('dest'));
                                    destPos = findstring(lqn.names,destName);
                                    destID = lqn.names{destPos};
                                    callPos = findstring(lqn.callnames,[actID,'=>',destID]);
                                    callResult = callElement.getElementsByTagName('result-call');
                                    wRes = str2double(callResult.item(0).getAttribute('waiting'));
                                    Avg.Edges.Waiting(callPos) = wRes;
                                end
                                %asynch-call
                                asynchCalls = actElement.getElementsByTagName('asynch-call');
                                for m = 0:asynchCalls.getLength()-1
                                    callElement = asynchCalls.item(m);
                                    destName = char(callElement.getAttribute('dest'));
                                    destPos = findstring(lqn.names,destName);
                                    destID = lqn.names{destPos};
                                    callPos = findstring(lqn.callnames,[actID,'->',destID]);
                                    callResult = callElement.getElementsByTagName('result-call');
                                    wRes = str2double(callResult.item(0).getAttribute('waiting'));
                                    Avg.Edges.Waiting(callPos) = wRes;
                                end
                            end
                        end
                    end
                end
            end
            
            self.result.RawAvg = Avg;
            self.result.Avg.ProcUtil = Avg.Nodes.ProcUtilization(:);
            self.result.Avg.SvcT = Avg.Nodes.Phase1ServiceTime(:);
            self.result.Avg.Tput = Avg.Nodes.Throughput(:);
            self.result.Avg.Util =  Avg.Nodes.Utilization(:);
            self.result.Avg.RespT = NaN*Avg.Nodes.ProcWaiting(:);
            self.result.Avg.SCV = Avg.Nodes.Phase1ServiceSCV(:);
            self.result.Avg.QLen = NaN*Avg.Nodes.ProcWaiting(:);
            result = self.result;
        end
 %------------------------------------------------------------------------------------------------       
        function [result, iterations, exceedProb] = parseXMLResultsExceedProb(self, filename)
            % [RESULT, ITERATIONS] = PARSEXMLRESULTS(FILENAME)
            
            import javax.xml.parsers.*;
            import org.w3c.dom.*;
            import java.io.*;
            
            exceedProb = [];
            
            lqn = self.getStruct;
            numOfNodes = lqn.nidx;
            numOfCalls = lqn.ncalls;
            Avg.Nodes.Utilization = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase1Utilization = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase2Utilization = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase1ServiceTime = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase1ServiceSCV = NaN*ones(numOfNodes,1);
            Avg.Nodes.Phase2ServiceTime = NaN*ones(numOfNodes,1);
            Avg.Nodes.Throughput = NaN*ones(numOfNodes,1);
            Avg.Nodes.ProcWaiting = NaN*ones(numOfNodes,1);
            Avg.Nodes.ProcUtilization = NaN*ones(numOfNodes,1);
            Avg.Edges.Waiting = NaN*ones(numOfCalls,1);
            
            % init Java XML parser and load file
            dbFactory = DocumentBuilderFactory.newInstance();
            dBuilder = dbFactory.newDocumentBuilder();
            
            [fpath,fname,~] = fileparts(filename);
            resultFilename = [fpath,filesep,fname,'.lqxo'];
            if self.options.verbose
                line_printf('\nParsing LQNS result file: %s',resultFilename);
            end
            if self.options.keep
                line_printf('\nLQNS result file available at: %s',resultFilename);
            end
            
            doc = dBuilder.parse(resultFilename);
            doc.getDocumentElement().normalize();
            
            %solver-params
            solverParams = doc.getElementsByTagName('solver-params');
            for i = 0:solverParams.getLength()-1
                solverParam = solverParams.item(i);
                result = solverParam.getElementsByTagName('result-general');
                iterations = str2double(result.item(0).getAttribute('iterations'));
            end
            
            procList = doc.getElementsByTagName('processor');
            for i = 0:procList.getLength()-1
                %Element - Host
                procElement = procList.item(i);
                procName = char(procElement.getAttribute('name'));
                procPos = findstring(lqn.names,procName);
                procResult = procElement.getElementsByTagName('result-processor');
                uRes = str2double(procResult.item(0).getAttribute('utilization'));
                Avg.Nodes.ProcUtilization(procPos) = uRes;
                
                taskList = procElement.getElementsByTagName('task');
                for j = 0:taskList.getLength()-1
                    %Element - Task
                    taskElement = taskList.item(j);
                    taskName = char(taskElement.getAttribute('name'));
                    taskPos = findstring(lqn.names,taskName);
                    taskResult = taskElement.getElementsByTagName('result-task');
                    uRes = str2double(taskResult.item(0).getAttribute('utilization'));
                    p1uRes = str2double(taskResult.item(0).getAttribute('phase1-utilization'));
                    p2uRes = str2double(taskResult.item(0).getAttribute('phase2-utilization'));
                    tRes = str2double(taskResult.item(0).getAttribute('throughput'));
                    puRes = str2double(taskResult.item(0).getAttribute('proc-utilization'));
                    Avg.Nodes.Utilization(taskPos) = uRes;
                    Avg.Nodes.Phase1Utilization(taskPos) = p1uRes;
                    Avg.Nodes.Phase2Utilization(taskPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                    Avg.Nodes.Throughput(taskPos) = tRes;
                    Avg.Nodes.ProcUtilization(taskPos) = puRes;
                    
                    entryList = taskElement.getElementsByTagName('entry');
                    for k = 0:entryList.getLength()-1
                        %Element - Entry
                        entryElement = entryList.item(k);
                        entryName = char(entryElement.getAttribute('name'));
                        entryPos = findstring(lqn.names,entryName);
                        entryResult = entryElement.getElementsByTagName('result-entry');
                        uRes = str2double(entryResult.item(0).getAttribute('utilization'));
                        p1uRes = str2double(entryResult.item(0).getAttribute('phase1-utilization'));
                        p1SCV = str2double(entryResult.item(0).getAttribute('squared-coeff-variation'));
                        p2uRes = str2double(entryResult.item(0).getAttribute('phase2-utilization'));
                        p1stRes = str2double(entryResult.item(0).getAttribute('phase1-service-time'));
                        p2stRes = str2double(entryResult.item(0).getAttribute('phase2-service-time'));
                        tRes = str2double(entryResult.item(0).getAttribute('throughput'));
                        puRes = str2double(entryResult.item(0).getAttribute('proc-utilization'));
                        Avg.Nodes.Utilization(entryPos) = uRes;
                        Avg.Nodes.Phase1Utilization(entryPos) = p1uRes;
                        Avg.Nodes.Phase2Utilization(entryPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                        Avg.Nodes.Phase1ServiceTime(entryPos) = p1stRes;
                        Avg.Nodes.Phase1ServiceSCV(entryPos) = p1SCV;
                        Avg.Nodes.Phase2ServiceTime(entryPos) = ifthenelse(isempty(p2stRes),NaN,p2stRes);
                        Avg.Nodes.Throughput(entryPos) = tRes;
                        Avg.Nodes.ProcUtilization(entryPos) = puRes;
                    end
                    
                    %task-activities
                    taskActsList = taskElement.getElementsByTagName('task-activities');
                    if taskActsList.getLength > 0
                        taskActsElement = taskActsList.item(0);
                        actList = taskActsElement.getElementsByTagName('activity');
                        for l = 0:actList.getLength()-1
                            %Element - Activity
                            actElement = actList.item(l);
                            if strcmp(char(actElement.getParentNode().getNodeName()),'task-activities')
                                actName = char(actElement.getAttribute('name'));
                                actPos = findstring(lqn.names,actName);
                                actResult = actElement.getElementsByTagName('result-activity');
                                uRes = str2double(actResult.item(0).getAttribute('utilization'));
                                stRes = str2double(actResult.item(0).getAttribute('service-time'));
                                tRes = str2double(actResult.item(0).getAttribute('throughput'));
                                pwRes = str2double(actResult.item(0).getAttribute('proc-waiting'));
                                puRes = str2double(actResult.item(0).getAttribute('proc-utilization'));
                                try
                                    exceedProb = [exceedProb;str2double(actResult.item(0).getAttribute('prob-exceed-max-service-time'))];
                                catch
                                    exceedProb;
                                end
                                Avg.Nodes.Utilization(actPos) = uRes;
                                Avg.Nodes.Phase1ServiceTime(actPos) = stRes;
                                Avg.Nodes.Throughput(actPos) = tRes;
                                Avg.Nodes.ProcWaiting(actPos) = pwRes;
                                Avg.Nodes.ProcUtilization(actPos) = puRes;
                                
                                actID = lqn.names{actPos};
                                %synch-call
                                synchCalls = actElement.getElementsByTagName('synch-call');
                                for m = 0:synchCalls.getLength()-1
                                    callElement = synchCalls.item(m);
                                    destName = char(callElement.getAttribute('dest'));
                                    destPos = findstring(lqn.names,destName);
                                    destID = lqn.names{destPos};
                                    callPos = findstring(lqn.callnames,[actID,'=>',destID]);
                                    callResult = callElement.getElementsByTagName('result-call');
                                    wRes = str2double(callResult.item(0).getAttribute('waiting'));
                                    Avg.Edges.Waiting(callPos) = wRes;
                                end
                                %asynch-call
                                asynchCalls = actElement.getElementsByTagName('asynch-call');
                                for m = 0:asynchCalls.getLength()-1
                                    callElement = asynchCalls.item(m);
                                    destName = char(callElement.getAttribute('dest'));
                                    destPos = findstring(lqn.names,destName);
                                    destID = lqn.names{destPos};
                                    callPos = findstring(lqn.callnames,[actID,'->',destID]);
                                    callResult = callElement.getElementsByTagName('result-call');
                                    wRes = str2double(callResult.item(0).getAttribute('waiting'));
                                    Avg.Edges.Waiting(callPos) = wRes;
                                end
                            end
                        end
                    end
                end
            end
            
            self.result.RawAvg = Avg;
            self.result.Avg.ProcUtil = Avg.Nodes.ProcUtilization(:);
            self.result.Avg.SvcT = Avg.Nodes.Phase1ServiceTime(:);
            self.result.Avg.Tput = Avg.Nodes.Throughput(:);
            self.result.Avg.Util =  Avg.Nodes.Utilization(:);
            self.result.Avg.RespT = NaN*Avg.Nodes.ProcWaiting(:);
            self.result.Avg.SCV = Avg.Nodes.Phase1ServiceSCV(:);
            self.result.Avg.QLen = NaN*Avg.Nodes.ProcWaiting(:);
            result = self.result;
            
            nanIndices = isnan(exceedProb);
            exceedProb(nanIndices) = [];
%             exceedProb = exceedProb(~nanIndices);
        end
%----------------------------------------------------------------------------------------------------        
    end
    
    methods %(Hidden)
        function [NodeAvgTable,CallAvgTable] = getRawAvgTables(self)
            % [QN,UN,RN,TN] = GETRAWAVGTABLES(SELF,~,~,~,~)
            
            runAnalyzer(self);
            
            lqn = self.getStruct;
            Node = label(lqn.names);
            O = length(Node);
            NodeType = label(O,1);
            for o = 1:O
                switch lqn.type(o)
                    case LayeredNetworkElement.PROCESSOR
                        NodeType(o,1) = label({'Processor'});
                    case LayeredNetworkElement.TASK
                        NodeType(o,1) = label({'Task'});
                    case LayeredNetworkElement.ENTRY
                        NodeType(o,1) = label({'Entry'});
                    case LayeredNetworkElement.ACTIVITY
                        NodeType(o,1) = label({'Activity'});
                    case LayeredNetworkElement.CALL
                        NodeType(o,1) = label({'Call'});
                end
            end
            Utilization = self.result.RawAvg.Nodes.Utilization;
            Phase1Utilization = self.result.RawAvg.Nodes.Phase1Utilization;
            Phase2Utilization = self.result.RawAvg.Nodes.Phase2Utilization;
            Phase1ServiceTime = self.result.RawAvg.Nodes.Phase1ServiceTime;
            Phase1ServiceSCV = self.result.RawAvg.Nodes.Phase1ServiceSCV;
            Phase2ServiceTime = self.result.RawAvg.Nodes.Phase2ServiceTime;
            Throughput = self.result.RawAvg.Nodes.Throughput;
            ProcWaiting = self.result.RawAvg.Nodes.ProcWaiting;
            ProcUtilization = self.result.RawAvg.Nodes.ProcUtilization;
            NodeAvgTable = Table(Node, NodeType, Utilization, Phase1Utilization,...
                Phase2Utilization, Phase1ServiceTime, Phase1ServiceSCV, Phase2ServiceTime, Throughput,...
                ProcWaiting, ProcUtilization);
            
            lqn=self.getStruct;
            %CallIndices = find(self.model.lqnGraph.Edges.Type>0);
            %CallIndices = lqn.callidx;
            %EndNodes = self.model.lqnGraph.Edges.EndNodes(CallIndices,:);
            %EndNodes = lqn.callpair(:,2);
            %SourceIndices = findnode(self.model.lqnGraph,EndNodes(:,1));
            %SourceNode = label(self.model.lqnGraph.Nodes.Node(SourceIndices));
            %TargetIndices = findnode(self.model.lqnGraph,EndNodes(:,2));
            %TargetNode = label(self.model.lqnGraph.Nodes.Node(TargetIndices));
            
            %CallTypeMap = [CallType.SYNC;CallType.ASYNC;CallType.FWD];
            %Type = label(CallTypeMap(self.model.lqnGraph.Edges.Type(CallIndices)));
            if lqn.ncalls == 0
                CallAvgTable = Table();
            else
                SourceNode = label({lqn.names{lqn.callpair(:,1)}})';
                TargetNode = label({lqn.names{lqn.callpair(:,2)}})';
                Type = lqn.calltype;
                Waiting = self.result.RawAvg.Edges.Waiting(1:lqn.ncalls);
                CallAvgTable = Table(SourceNode, TargetNode, Type, Waiting);
            end
        end
%---------------------------------------------------------------------------------------------------        
        function [Phase1ServiceTime,Phase1ServiceSCV,runningtime] = getRawTimes(self)
            % [QN,UN,RN,TN] = GETRAWAVGTABLES(SELF,~,~,~,~)
            
            runningtime = runAnalyzer(self);
            
            lqn = self.getStruct;
            Node = label(lqn.names);
            O = length(Node);
            NodeType = label(O,1);
            for o = 1:O
                switch lqn.type(o)
                    case LayeredNetworkElement.PROCESSOR
                        NodeType(o,1) = label({'Processor'});
                    case LayeredNetworkElement.TASK
                        NodeType(o,1) = label({'Task'});
                    case LayeredNetworkElement.ENTRY
                        NodeType(o,1) = label({'Entry'});
                    case LayeredNetworkElement.ACTIVITY
                        NodeType(o,1) = label({'Activity'});
                    case LayeredNetworkElement.CALL
                        NodeType(o,1) = label({'Call'});
                end
            end
            Phase1ServiceTime = self.result.RawAvg.Nodes.Phase1ServiceTime;
            Phase1ServiceSCV = self.result.RawAvg.Nodes.Phase1ServiceSCV;
        end
%----------------------------------------------------------------------------------------------------        
        function [Phase1ServiceTime,Phase1ServiceSCV,runningtime,exceedProb] = getRawTimes_thresholdTime(self,thresholdTime)
            % [QN,UN,RN,TN] = GETRAWAVGTABLES(SELF,~,~,~,~)
            
            [runningtime,exceedProb] = runAnalyzer_thresholdTime(self, thresholdTime);
            
            lqn = self.getStruct;
            Node = label(lqn.names);
            O = length(Node);
            NodeType = label(O,1);
            for o = 1:O
                switch lqn.type(o)
                    case LayeredNetworkElement.PROCESSOR
                        NodeType(o,1) = label({'Processor'});
                    case LayeredNetworkElement.TASK
                        NodeType(o,1) = label({'Task'});
                    case LayeredNetworkElement.ENTRY
                        NodeType(o,1) = label({'Entry'});
                    case LayeredNetworkElement.ACTIVITY
                        NodeType(o,1) = label({'Activity'});
                    case LayeredNetworkElement.CALL
                        NodeType(o,1) = label({'Call'});
                end
            end
            Phase1ServiceTime = self.result.RawAvg.Nodes.Phase1ServiceTime;
            Phase1ServiceSCV = self.result.RawAvg.Nodes.Phase1ServiceSCV;
        end
%----------------------------------------------------------------------------------------------------
    end
    
    methods (Static)
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Queue',...
                'Coxian','Erlang','Exponential','HyperExp',...
                'Buffer','Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_PS','SchedStrategy_FCFS','ClosedClass'});
            bool = true;
            for e=1:model.getNumberOfLayers()
                bool = bool && SolverFeatureSet.supports(featSupported, featUsed{e});
            end
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('LQNS');
        end
        
        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            
            bool = true;
            if ispc % windows
                [~,ret] = dos('lqns -V -H');
                if containsstr(ret,'not recognized')
                    bool = false;
                end
                if containsstr(ret,'Version 5') || containsstr(ret,'Version 4') ...
                        || containsstr(ret,'Version 3') || containsstr(ret,'Version 2') ...
                        || containsstr(ret,'Version 1')
                    line_warning(mfilename,'Unsupported LQNS version. LINE requires Version 6.0 or greater.');
                    bool = true;
                end
            else % linux
                [~,ret] = unix('lqns -V -H');
                if containsstr(ret,'command not found')
                    bool = false;
                end
                if containsstr(ret,'Version 5') || containsstr(ret,'Version 4') ...
                        || containsstr(ret,'Version 3') || containsstr(ret,'Version 2') ...
                        || containsstr(ret,'Version 1')
                    line_warning(mfilename,'Unsupported LQNS version. LINE requires Version 6.0 or greater.');
                    bool = true;
                end
            end
        end
        
        function plot(model)
            lastwd = pwd;
            cwd = [lineRootFolder,filesep,'workspace',filesep,'lqns'];
            cd(cwd)
            fname = tempname([lineRootFolder,filesep,'workspace',filesep,'lqns']);
            model.writeXML(fname);
            [d,f]=fileparts(fname);
            system(sprintf('lqn2ps %s',f));
            fprintf(1,'Postscript file saved in: %s\n',[fname,'.ps'])
        end
    end

end

