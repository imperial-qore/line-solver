function [result, iterations] = parseXMLResults(self, filename)
% [RESULT, ITERATIONS] = PARSEXMLRESULTS(FILENAME)

import javax.xml.parsers.*; %#ok<JAPIEXT629>
import org.w3c.dom.*;        %#ok<JAPIEXT161>
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
    if (procNode.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE) %#ok<JAPIEXT161>
        procElement = procNode;
        procName = char(procElement.getAttribute('name'));
        procPos = findlqnelem(lqn,procName,LayeredNetworkElement.HOST);
        procResult = procElement.getElementsByTagName('result-processor');
        uRes = str2double(procResult.item(0).getAttribute('utilization'));
        if procPos > 0
            Avg.Nodes.ProcUtilization(procPos) = uRes;
        end

        taskList = procElement.getElementsByTagName('task');
        for j = 0:taskList.getLength()-1
            %Element - Task
            taskElement = taskList.item(j);
            taskName = char(taskElement.getAttribute('name'));
            taskPos = findlqnelem(lqn,taskName,LayeredNetworkElement.TASK);
            taskResult = taskElement.getElementsByTagName('result-task');
            uRes = str2double(taskResult.item(0).getAttribute('utilization'));
            p1uRes = str2double(taskResult.item(0).getAttribute('phase1-utilization'));
            p2uRes = str2double(taskResult.item(0).getAttribute('phase2-utilization'));
            tRes = str2double(taskResult.item(0).getAttribute('throughput'));
            puRes = str2double(taskResult.item(0).getAttribute('proc-utilization'));
            if taskPos > 0
                Avg.Nodes.Utilization(taskPos) = uRes;
                Avg.Nodes.Phase1Utilization(taskPos) = p1uRes;
                Avg.Nodes.Phase2Utilization(taskPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                Avg.Nodes.Throughput(taskPos) = tRes;
                Avg.Nodes.ProcUtilization(taskPos) = puRes;
            end

            entryList = taskElement.getElementsByTagName('entry');
            for k = 0:entryList.getLength()-1
                %Element - Entry
                entryElement = entryList.item(k);
                entryName = char(entryElement.getAttribute('name'));
                entryPos = findlqnelem(lqn,entryName,LayeredNetworkElement.ENTRY);
                entryResult = entryElement.getElementsByTagName('result-entry');
                uRes = str2double(entryResult.item(0).getAttribute('utilization'));
                p1uRes = str2double(entryResult.item(0).getAttribute('phase1-utilization'));
                p2uRes = str2double(entryResult.item(0).getAttribute('phase2-utilization'));
                p1stRes = str2double(entryResult.item(0).getAttribute('phase1-service-time'));
                p2stRes = str2double(entryResult.item(0).getAttribute('phase2-service-time'));
                tRes = str2double(entryResult.item(0).getAttribute('throughput'));
                puRes = str2double(entryResult.item(0).getAttribute('proc-utilization'));
                if entryPos > 0
                    Avg.Nodes.Utilization(entryPos) = uRes;
                    Avg.Nodes.Phase1Utilization(entryPos) = p1uRes;
                    Avg.Nodes.Phase2Utilization(entryPos) = ifthenelse(isempty(p2uRes),NaN,p2uRes);
                    Avg.Nodes.Phase1ServiceTime(entryPos) = p1stRes;
                    Avg.Nodes.Phase2ServiceTime(entryPos) = ifthenelse(isempty(p2stRes),NaN,p2stRes);
                    Avg.Nodes.Throughput(entryPos) = tRes;
                    Avg.Nodes.ProcUtilization(entryPos) = puRes;
                end

                %entry-phase-activities (PH1PH2 format): fill the phase
                %activity rows, which are otherwise left NaN (JLINE parity)
                epaList = entryElement.getElementsByTagName('entry-phase-activities');
                if epaList.getLength > 0
                    epaElement = epaList.item(0);
                    phActList = epaElement.getElementsByTagName('activity');
                    for l = 0:phActList.getLength()-1
                        actElement = phActList.item(l);
                        actName = char(actElement.getAttribute('name'));
                        actPos = findlqnelem(lqn,actName,LayeredNetworkElement.ACTIVITY);
                        if actPos <= 0
                            continue;
                        end
                        actResult = actElement.getElementsByTagName('result-activity');
                        if actResult.getLength == 0
                            continue;
                        end
                        uResAct = str2double(actResult.item(0).getAttribute('utilization'));
                        stRes = str2double(actResult.item(0).getAttribute('service-time'));
                        tResAct = str2double(actResult.item(0).getAttribute('throughput'));
                        pwRes = str2double(actResult.item(0).getAttribute('proc-waiting'));
                        puResAct = str2double(actResult.item(0).getAttribute('proc-utilization'));
                        hdRes = str2double(actElement.getAttribute('host-demand-mean'));
                        Avg.Nodes.Utilization(actPos) = uResAct;
                        Avg.Nodes.Phase1ServiceTime(actPos) = stRes;
                        Avg.Nodes.ProcWaiting(actPos) = pwRes;
                        if ~isnan(tResAct)
                            Avg.Nodes.Throughput(actPos) = tResAct;
                        else
                            % LQNS omits throughput from entry-phase result-activity
                            % elements; each phase executes once per entry invocation
                            Avg.Nodes.Throughput(actPos) = tRes;
                        end
                        if ~isnan(puResAct)
                            Avg.Nodes.ProcUtilization(actPos) = puResAct;
                        elseif ~isnan(tRes) && ~isnan(hdRes)
                            % LQNS omits proc-utilization from entry-phase
                            % result-activity elements; the per-phase value is
                            % the entry throughput times the phase host demand
                            Avg.Nodes.ProcUtilization(actPos) = tRes * hdRes;
                        end
                        actID = lqn.names{actPos};
                        %synch-call
                        synchCalls = actElement.getElementsByTagName('synch-call');
                        for m = 0:synchCalls.getLength()-1
                            callElement = synchCalls.item(m);
                            destName = char(callElement.getAttribute('dest'));
                            destPos = findlqnelem(lqn,destName,LayeredNetworkElement.ENTRY);
                            if destPos <= 0
                                continue;
                            end
                            destID = lqn.names{destPos};
                            callPos = findstring(lqn.callnames,[actID,'=>',destID]);
                            callResult = callElement.getElementsByTagName('result-call');
                            if callPos > 0 && callResult.getLength > 0
                                wRes = str2double(callResult.item(0).getAttribute('waiting'));
                                Avg.Edges.Waiting(callPos) = wRes;
                            end
                        end
                        %asynch-call
                        asynchCalls = actElement.getElementsByTagName('asynch-call');
                        for m = 0:asynchCalls.getLength()-1
                            callElement = asynchCalls.item(m);
                            destName = char(callElement.getAttribute('dest'));
                            destPos = findlqnelem(lqn,destName,LayeredNetworkElement.ENTRY);
                            if destPos <= 0
                                continue;
                            end
                            destID = lqn.names{destPos};
                            callPos = findstring(lqn.callnames,[actID,'->',destID]);
                            callResult = callElement.getElementsByTagName('result-call');
                            if callPos > 0 && callResult.getLength > 0
                                wRes = str2double(callResult.item(0).getAttribute('waiting'));
                                Avg.Edges.Waiting(callPos) = wRes;
                            end
                        end
                    end
                end
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
                        actPos = findlqnelem(lqn,actName,LayeredNetworkElement.ACTIVITY);
                        if actPos <= 0
                            continue;
                        end
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
                            destPos = findlqnelem(lqn,destName,LayeredNetworkElement.ENTRY);
                            if destPos <= 0
                                continue;
                            end
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
                            destPos = findlqnelem(lqn,destName,LayeredNetworkElement.ENTRY);
                            if destPos <= 0
                                continue;
                            end
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

% Processor utilization of an entry, aggregated from its activity graph.
% lqns credits host work to whichever level carries the host demand: in the
% activity-graph form -- the only form writeXML emits -- an entry declares
% none, so lqns reports result-entry proc-utilization as a literal 0 and the
% work sits on the result-activity rows. The entry value is then the sum over
% the activities reachable from the entry within its own task, which is what
% lqn.actsof holds. In PH1PH2 form the same sum runs over the phase activities
% and reproduces the value lqns reports there, so no form test is needed. An
% entry with no activities, or any activity lqns left unreported, keeps the
% raw attribute rather than a partial sum.
for eoff = 1:lqn.nentries
    eidx = lqn.eshift + eoff;
    acts = lqn.actsof{eidx};
    if isempty(acts)
        continue
    end
    puActs = Avg.Nodes.ProcUtilization(acts);
    if ~any(isnan(puActs))
        Avg.Nodes.ProcUtilization(eidx) = sum(puActs);
    end
end

% Phase-1 service time of an entry lqns never invoked.
% lqns omits phase1-service-time from result-entry exactly when the entry's
% throughput is zero: nothing was served, so there is no per-invocation mean to
% report. LINE then carried a NaN where the table says an entry HAS a response
% time and every other solver reports one, breaking the NaN mask -- see
% _kb/06-solver-catalog.md. The value is taken from the activity rows, and ONLY
% where they are unanimous: if every activity reachable from the entry reports a
% zero service time then every aggregation law agrees on zero -- the serial sum,
% the branch-weighted mean of an OrFork, the order statistic of an AndFork -- so
% the derivation does not depend on which one applies.
% It is deliberately NOT generalised the way ProcUtilization is above.
% Utilizations add over an activity graph; response times do not. Measured over
% the example corpus, sum(actsof) reproduces phase1-service-time on serial
% chains only and misses it wherever the graph branches (lqn_workflows `Entry`:
% 12.5667 reported against 8.5667 summed, lqn_fork_open_arrival `SE`: 0.841667
% against 1.0), so a summed fallback would answer with a number lqns contradicts.
% An entry whose activities are unreported, absent, or not all zero keeps NaN.
for eoff = 1:lqn.nentries
    eidx = lqn.eshift + eoff;
    if ~isnan(Avg.Nodes.Phase1ServiceTime(eidx))
        continue
    end
    acts = lqn.actsof{eidx};
    if isempty(acts)
        continue
    end
    stActs = Avg.Nodes.Phase1ServiceTime(acts);
    if ~any(isnan(stActs)) && all(stActs == 0)
        Avg.Nodes.Phase1ServiceTime(eidx) = 0;
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
