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
    line_printf('\nParsing LQNS result file: %s\n',resultFilename);
    if self.options.keep
        %line_printf('\nLQNS result file available at: %s',resultFilename);
    end
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
    procNode = procList.item(i);
    if (procNode.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE)
        procElement = procNode;
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
                p2uRes = str2double(entryResult.item(0).getAttribute('phase2-utilization'));
                p1stRes = str2double(entryResult.item(0).getAttribute('phase1-service-time'));
                p2stRes = str2double(entryResult.item(0).getAttribute('phase2-service-time'));
                tRes = str2double(entryResult.item(0).getAttribute('throughput'));
                puRes = str2double(entryResult.item(0).getAttribute('proc-utilization'));
                Avg.Nodes.Utilization(entryPos) = uRes;
                Avg.Nodes.Phase1Utilization(entryPos) = p1uRes;
                Avg.Nodes.Phase2Utilization(entryPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                Avg.Nodes.Phase1ServiceTime(entryPos) = p1stRes;
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
end

self.result.RawAvg = Avg;
self.result.Avg.ProcUtil = Avg.Nodes.ProcUtilization(:);
self.result.Avg.SvcT = Avg.Nodes.Phase1ServiceTime(:);
self.result.Avg.Tput = Avg.Nodes.Throughput(:);
self.result.Avg.Util =  Avg.Nodes.Utilization(:);
self.result.Avg.RespT = NaN*Avg.Nodes.ProcWaiting(:);
self.result.Avg.QLen = NaN*Avg.Nodes.ProcWaiting(:);
result = self.result;
end
