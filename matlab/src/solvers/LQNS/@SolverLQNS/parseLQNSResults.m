function [Gres, iterations] =parseLQNSResults(G, filename, filename_sim, verbose)
% [GRES, ITERATIONS] =PARSELQNSRESULTS(G, FILENAME, FILENAME_SIM, VERBOSE)

% Copyright (c) 2012-2021, Imperial College London
% All rights reserved.

[Gres, iterations] = aux_results(G, filename, verbose);
Gres.Nodes = movevars(Gres.Nodes,'Wlqns','Before','U');
Gres.Nodes = movevars(Gres.Nodes,'Ulqns','Before','T');
Gres.Nodes = movevars(Gres.Nodes,'Tlqns','Before','Q');
Gres.Nodes.RespT(Gres.Nodes.RespT<1e-4)=0;
Gres.Nodes.Util(Gres.Nodes.Util<1e-4)=0;
Gres.Nodes.Tput(Gres.Nodes.Tput<1e-4)=0;
Gres.Nodes.QLen(Gres.Nodes.QLen<1e-4)=0;
Gres.Nodes.Name = [];

Gsim = aux_results(G, filename_sim, false);
Gsim.Nodes = movevars(Gsim.Nodes,'Wlqns','Before','U');
Gsim.Nodes = movevars(Gsim.Nodes,'Ulqns','Before','T');
Gsim.Nodes = movevars(Gsim.Nodes,'Tlqns','Before','Q');
Gsim.Nodes.RespT(Gsim.Nodes.RespT<1e-4)=0;
Gsim.Nodes.Util(Gsim.Nodes.Util<1e-4)=0;
Gsim.Nodes.Tput(Gsim.Nodes.Tput<1e-4)=0;
Gsim.Nodes.QLen(Gsim.Nodes.QLen<1e-4)=0;
Gsim.Nodes.Name = [];
Gres.Nodes.RespTsim = Gsim.Nodes.RespTlqns;
Gres.Nodes.Utilsim = Gsim.Nodes.Utillqns;
Gres.Nodes.Tputsim = Gsim.Nodes.Tputlqns;
Gres.Nodes.QLensim = Gsim.Nodes.QLenlqns;
Gres.Nodes = movevars(Gres.Nodes,'Wsim','Before','U');
Gres.Nodes = movevars(Gres.Nodes,'Usim','Before','T');
Gres.Nodes = movevars(Gres.Nodes,'Tsim','Before','Q');
end

function [Gres, iterations] = aux_results(G, filename, verbose)
% [GRES, ITERATIONS] = AUX_RESULTS(G, FILENAME, VERBOSE)

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;

import LayeredNetwork.*;

Gres = G;
Gres.Nodes.Type = [];
Gres.Nodes.Proc = [];
Gres.Nodes.Task = [];
Gres.Nodes.Entry = [];
Gres.Nodes.D = [];
Gres.Nodes.Object = [];
%Gres.Nodes.Node = [];
Gres.Nodes.MaxJobs = [];
Gres.Nodes.Mult = [];
Gres.Nodes.RespTlqns(:) = 0;
Gres.Nodes.Utillqns(:) = 0;
Gres.Nodes.Tputlqns(:) = 0;
Gres.Nodes.QLenlqns(:) = 0;

% LQN
myLN = LayeredNetwork(filename);

if nargin == 1
    verbose = 0;
end

% init Java XML parser and load file
dbFactory = DocumentBuilderFactory.newInstance();
dBuilder = dbFactory.newDocumentBuilder();
try
    doc = dBuilder.parse(filename);
catch exception %java.io.FileNotFoundException
    if ~exist(filename, 'file')
        line_printf(['Error: Input XML file ', filename, ' not found']);
        hosts = [];
        tasks = [];
        entries = [];
        requesters = [];
        %providers = [];
        myLN.hosts = hosts;
        %myLN.tasks = tasks;
        %myLN.entries = entries;
        
        return;
    else
        rethrow(exception);
    end
end

doc.getDocumentElement().normalize();
if verbose > 0
    line_printf(['Parsing LQN file: ', filename] );
    line_printf(['Root element :', char(doc.getDocumentElement().getNodeName()) ] );
end

%NodeList
solverPara = doc.getElementsByTagName('solver-params');
for i = 0:solverPara.getLength()-1
    procNode = solverPara.item(i);
    result = procNode.getElementsByTagName('result-general');
    iterations = str2num(result.item(0).getAttribute('iterations'));
end

procList = doc.getElementsByTagName('processor');
hosts = [];
%providers = cell(0); % list of entries that provide services - Entry, Task, Proc
requesters = cell(0); % list of activities that request services - Act, Task, Proc
tasks = cell(0); %list of tasks - Task, task ID, Proc, ProcID - Row Index as task ID
entries = cell(0); %list of entries - Entry, Task ID
taskID = 1;
physical = cell(0); %list of actual hosts, those thata receive demand for resources
%demand is always indicated in an entry activity
procID = 1;

clients = []; % list of tasks that act as pure clients (think time)
for i = 0:procList.getLength()-1
    %Node - Host
    procNode = procList.item(i);
    
    if procNode.getNodeType() == Node.ELEMENT_NODE
        
        %Element
        procElement = procNode;
        name = char(procElement.getAttribute('name'));
        result = procNode.getElementsByTagName('result-processor');
        Gres.Nodes.Tputlqns(findstring(G.Nodes.Node,name)) = 0;
        utilizationRes = str2num(result.item(0).getAttribute('utilization'));
        Gres.Nodes.Utillqns(findstring(G.Nodes.Node,name)) = utilizationRes;
        
        taskList = procNode.getElementsByTagName('task');
        for j = 0:taskList.getLength()-1
            %Node - Task
            taskNode = taskList.item(j);
            if taskNode.getNodeType() == Node.ELEMENT_NODE
                %Element
                taskElement = taskNode;
                name = char(taskElement.getAttribute('name'));
                result = taskNode.getElementsByTagName('result-task');
                utilizationRes = str2num(result.item(0).getAttribute('proc-utilization'));
                Gres.Nodes.Utillqns(findstring(G.Nodes.Node,name)) = utilizationRes;
                qlenRes = str2num(result.item(0).getAttribute('utilization'));
                Gres.Nodes.QLenlqns(findstring(G.Nodes.Node,name)) = qlenRes;
                tputRes = str2num(result.item(0).getAttribute('throughput'));
                Gres.Nodes.Tputlqns(findstring(G.Nodes.Node,name)) = tputRes;                  
                entryList = taskNode.getElementsByTagName('entry');
                for k = 0:entryList.getLength()-1
                    %Node - Task
                    entryNode = entryList.item(k);
                    if entryNode.getNodeType() == Node.ELEMENT_NODE
                        %Element
                        entryElement = entryNode;
                        name = char(entryElement.getAttribute('name'));
                        result = entryNode.getElementsByTagName('result-entry');
                        utilizationRes = str2num(result.item(0).getAttribute('proc-utilization'));
                        Gres.Nodes.Utillqns(findstring(G.Nodes.Node,name)) = utilizationRes*NaN; % ignore as this different in LINE
                        qlenRes = str2num(result.item(0).getAttribute('utilization'));
                        Gres.Nodes.QLenlqns(findstring(G.Nodes.Node,name)) = qlenRes;
                        tputRes = str2num(result.item(0).getAttribute('throughput'));
                        Gres.Nodes.Tputlqns(findstring(G.Nodes.Node,name)) = tputRes;
                        rtRes = str2num(result.item(0).getAttribute('phase1-service-time'));
                        if ~isempty(rtRes)
                            Gres.Nodes.RespTlqns(findstring(G.Nodes.Node,name)) = rtRes;
                        end
                    end
                end
                
                %% task-activities
                if taskElement.getElementsByTagName('task-activities').getLength > 0
                    %actNames = cell(0); iterActNames = 1;
                    %actCalls = cell(0);
                    actList = taskElement.getElementsByTagName('task-activities').item(0).getElementsByTagName('activity');
                    for l = 0:actList.getLength()-1
                        %Node - Task
                        actNode = actList.item(l);
                        if actNode.getNodeType() == Node.ELEMENT_NODE && strcmp(char(actNode.getParentNode().getNodeName()),'task-activities')
                            %Element
                            actElement = actNode;
                            name = char(actElement.getAttribute('name'));
                            
                            result = actNode.getElementsByTagName('result-activity');
                            rtRes = str2num(result.item(0).getAttribute('service-time'));
                            Gres.Nodes.RespTlqns(findstring(G.Nodes.Node,name)) = rtRes;
                            utilizationRes = str2num(result.item(0).getAttribute('proc-utilization'));
                            Gres.Nodes.Utillqns(findstring(G.Nodes.Node,name)) = utilizationRes;
                            qlenRes = str2num(result.item(0).getAttribute('utilization'));
                            Gres.Nodes.QLenlqns(findstring(G.Nodes.Node,name)) = qlenRes;
                            tputRes = str2num(result.item(0).getAttribute('throughput'));
                            Gres.Nodes.Tputlqns(findstring(G.Nodes.Node,name)) = tputRes;
                            
                            %                             synchCalls = actElement.getElementsByTagName('synch-call');
                            %                             asynchCalls = actElement.getElementsByTagName('asynch-call');
                            %                             %add synch calls if any
                            %                             if synchCalls.getLength() > 0
                            %                                 for m = 0:synchCalls.getLength()-1
                            %                                     callElement = synchCalls.item(m);
                            %                                     dest = char(callElement.getAttribute('dest'));
                            %                                     mean = str2double(char(callElement.getAttribute('calls-mean')));
                            %                                     tempAct = tempAct.synchCall(dest,mean);
                            %                                     actCalls{iterActNames,1} = dest;
                            %                                     requesters{size(requesters,1)+1,1} = tempAct.name;
                            %                                     requesters{size(requesters,1),2} = taskID;
                            %                                     requesters{size(requesters,1),3} = tempProc.name;
                            %                                     requesters{size(requesters,1),4} = dest;
                            %                                     requesters{size(requesters,1),5} = procID;
                            %                                     %requesters:
                            %                                     % activity - task - processor - dest (entry) - procID
                            %                                 end
                            %                                 %else
                            %                                 %    actCalls{iterActNames,1} = [];
                            %                                 %end
                            %                                 %iterActNames = iterActNames + 1;
                            %                                 %add asynch calls if any
                            %                             elseif asynchCalls.getLength() > 0
                            %                                 for m = 0:asynchCalls.getLength()-1
                            %                                     callElement = asynchCalls.item(m);
                            %                                     dest = char(callElement.getAttribute('dest'));
                            %                                     mean = str2double(char(callElement.getAttribute('calls-mean')));
                            %                                     tempAct = tempAct.asynchCall(dest,mean);
                            %                                     actCalls{iterActNames,1} = dest;
                            %                                     requesters{size(requesters,1)+1,1} = tempAct.name;
                            %                                     requesters{size(requesters,1),2} = taskID;
                            %                                     requesters{size(requesters,1),3} = tempProc.name;
                            %                                     requesters{size(requesters,1),4} = dest;
                            %                                     requesters{size(requesters,1),5} = procID;
                            %                                 end
                            %                             else
                            %                                 actCalls{iterActNames,1} = [];
                            %                             end
                            %iterActNames = iterActNames + 1;
                        end
                    end
                end
            end
        end
    end
end
for edge=1:height(G.Edges)
    if G.Edges.Type(edge) == 1 % add contribution of sync-calls
        syncSource = G.Edges.EndNodes{edge,1};
        aidx = findstring(Gres.Nodes.Name,syncSource);
        if G.Edges.Weight(edge) >= 1
            Gres.Nodes.RespT(aidx) = Gres.Nodes.RespT(aidx) +  G.Edges.RespT(edge) * G.Edges.Weight(edge);
        else
            Gres.Nodes.RespT(aidx) = Gres.Nodes.RespT(aidx) +  G.Edges.RespT(edge);
        end
    end
end
end

