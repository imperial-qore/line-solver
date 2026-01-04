function myLN = parseXML(filename, verbose)
% MYLN = PARSEXML(FILENAME, VERBOSE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import java.io.File;

import LayeredNetwork.*;

% LQN
myLN = LayeredNetwork(strrep(filename,'_','\_'));

if nargin<2%~exist('verbose','var')
    verbose = 0;
end

% init Java XML parser and load file
dbFactory = DocumentBuilderFactory.newInstance();
dBuilder = dbFactory.newDocumentBuilder();

fid=fopen(filename,'r');
if fid==-1
    line_error(mfilename,'File cannot be found. Verify the current directory and the specified filename.');
else
    fclose(fid);
end

if isempty(fileparts(filename))
    doc = dBuilder.parse(which(filename));
else
    doc = dBuilder.parse(filename);
end
doc.getDocumentElement().normalize();
if verbose > 0
    line_printf(['Parsing LQN file: ',filename]);
    line_printf(['Root element :',char(doc.getDocumentElement().getNodeName())]);
end

hosts = cell(0); %list of hosts - Proc
tasks = cell(0); %list of tasks - Task, ProcID
entries = cell(0); %list of entries - Entry, TaskID, ProcID
activities = cell(0); %list of activities - Act, TaskID, ProcID
procID = 1;
taskID = 1;
entryID = 1;
actID = 1;
procObj = cell(0);
taskObj = cell(0);
entryObj = cell(0);
actObj = cell(0);

procList = doc.getElementsByTagName('processor');
for i = 0:procList.getLength()-1
    %Element - Host
    procElement = procList.item(i);
    name = char(procElement.getAttribute('name'));
    scheduling = char(procElement.getAttribute('scheduling'));
    multiplicity = str2double(char(procElement.getAttribute('multiplicity')));
    replication = str2double(char(procElement.getAttribute('replication')));
    
    if isnan(replication)
        replication=1;
    end
    if strcmp(scheduling, 'inf')
        if isfinite(multiplicity)
            line_warning(mfilename,'A finite multiplicity is specified for a host processor with INF scheduling. Remove it or set it to "inf".\n');
        end
        multiplicity = Inf;
    elseif isnan(multiplicity)
        multiplicity = 1;
    end
    quantum = str2double(char(procElement.getAttribute('quantum')));
    if isnan(quantum)
        quantum = 0.001;
    end
    speedFactor = str2double(char(procElement.getAttribute('speed-factor')));
    if isnan(speedFactor)
        speedFactor = 1.0;
    end
    newProc = Processor(myLN, name, multiplicity, SchedStrategy.fromText(scheduling), quantum, speedFactor);
    newProc.setReplication(replication);
    procObj{end+1,1} = newProc;
    
    taskList = procElement.getElementsByTagName('task');
    for j = 0:taskList.getLength()-1
        %Element - Task
        taskElement = taskList.item(j);
        name = char(taskElement.getAttribute('name'));
        scheduling = char(taskElement.getAttribute('scheduling'));
        replication = str2double(char(taskElement.getAttribute('replication')));
        if isnan(replication)
            replication=1;
        end
        
        multiplicity = str2double(char(taskElement.getAttribute('multiplicity')));
        if strcmp(scheduling, 'inf')
            if isfinite(multiplicity) 
                line_warning(mfilename,'A finite multiplicity is specified for a task with inf scheduling. Remove it or set it to inf.\n');
            end
            multiplicity = Inf;
        elseif isnan(multiplicity)
            multiplicity = 1;
        end
        thinkTimeMean = str2double(char(taskElement.getAttribute('think-time')));
        if isnan(thinkTimeMean)
            thinkTimeMean = 0.0;
        end
        if thinkTimeMean <= 0.0
            thinkTime = Immediate.getInstance();
        else
            thinkTime = Exp.fitMean(thinkTimeMean);
        end
        newTask = Task(myLN, name, multiplicity, SchedStrategy.fromText(scheduling), thinkTime);
        newTask.setReplication(replication);

        % Parse priority attribute if present
        priorityStr = char(taskElement.getAttribute('priority'));
        if ~isempty(priorityStr)
            priority = str2double(priorityStr);
            if ~isnan(priority)
                newTask.setPriority(priority);
            end
        end

        % Parse fan-in element if present (used for replication load distribution)
        fanInList = taskElement.getElementsByTagName('fan-in');
        if fanInList.getLength() > 0
            fanInElement = fanInList.item(0);
            source = char(fanInElement.getAttribute('source'));
            valueStr = char(fanInElement.getAttribute('value'));
            if ~isempty(source) && ~isempty(valueStr)
                value = str2double(valueStr);
                if ~isnan(value)
                    newTask.setFanIn(source, value);
                end
            end
        end

        taskObj{end+1,1} = newTask;

        entryList = taskElement.getElementsByTagName('entry');
        for k = 0:entryList.getLength()-1
            %Element - Entry
            entryElement = entryList.item(k);
            name = char(entryElement.getAttribute('name'));
            newEntry = Entry(myLN, name);
            openArrivalRate = str2double(char(entryElement.getAttribute('open-arrival-rate')));
            if ~isnan(openArrivalRate)
                newEntry.openArrivalRate = openArrivalRate;
            end

            % Parse entry type attribute
            eType = char(entryElement.getAttribute('type'));
            if ~isempty(eType)
                newEntry.setType(eType);
            end

            entryObj{end+1,1} = newEntry;

            % Parse forwarding calls
            forwardingList = entryElement.getElementsByTagName('forwarding');
            for fw = 0:forwardingList.getLength()-1
                fwdElement = forwardingList.item(fw);
                destName = char(fwdElement.getAttribute('dest'));
                probStr = char(fwdElement.getAttribute('prob'));
                if isempty(probStr)
                    prob = 1.0;
                else
                    prob = str2double(probStr);
                end
                newEntry.forward(destName, prob);
            end

            %entry-phase-activities
            entryPhaseActsList = entryElement.getElementsByTagName('entry-phase-activities');
            if entryPhaseActsList.getLength > 0
                entryPhaseActsElement = entryPhaseActsList.item(0);
                actList = entryPhaseActsElement.getElementsByTagName('activity');
                name = cell(actList.getLength(),1);
                for l = 0:actList.getLength()-1
                    %Element - Activity
                    actElement = actList.item(l);
                    phase = str2double(char(actElement.getAttribute('phase')));
                    name{phase} = char(actElement.getAttribute('name'));
                    hostDemandMean = str2double(char(actElement.getAttribute('host-demand-mean')));
                    hostDemandSCV = str2double(char(actElement.getAttribute('host-demand-cvsq')));
                    if isnan(hostDemandSCV)
                        hostDemandSCV = 1.0;
                    end
                    if hostDemandMean <= 0.0
                        hostDemand = Immediate.getInstance();
                    else
                        if hostDemandSCV <= 0.0
                            hostDemand = Det(hostDemandMean);
                        elseif hostDemandSCV == 1.0
                            hostDemand = Exp.fitMean(hostDemandMean);
                        else
                            hostDemand = APH.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                        end
                    end
                    if phase == 1
                        boundToEntry = newEntry.name;
                    else
                        boundToEntry = '';
                    end
                    callOrder = char(actElement.getAttribute('call-order'));
                    newAct = Activity(myLN, name{phase}, hostDemand, boundToEntry, callOrder);
                    newAct.phase = phase;  % Store the phase number

                    % Parse activity think-time
                    actThinkTimeMean = str2double(char(actElement.getAttribute('think-time')));
                    if ~isnan(actThinkTimeMean) && actThinkTimeMean > 0.0
                        newAct.setThinkTime(actThinkTimeMean);
                    end

                    actObj{end+1,1} = newAct;

                    %synch-call
                    synchCalls = actElement.getElementsByTagName('synch-call');
                    for m = 0:synchCalls.getLength()-1
                        callElement = synchCalls.item(m);
                        dest = char(callElement.getAttribute('dest'));
                        mean = str2double(char(callElement.getAttribute('calls-mean')));
                        newAct = newAct.synchCall(dest,mean);
                    end
                    
                    %asynch-call
                    asynchCalls = actElement.getElementsByTagName('asynch-call');
                    for m = 0:asynchCalls.getLength()-1
                        callElement = asynchCalls.item(m);
                        dest = char(callElement.getAttribute('dest'));
                        mean = str2double(char(callElement.getAttribute('calls-mean')));
                        newAct = newAct.asynchCall(dest,mean);
                    end
                    
                    activities{end+1,1} = newAct.name;
                    activities{end,2} = taskID;
                    activities{end,3} = procID;
                    newTask = newTask.addActivity(newAct);
                    newAct.parent = newTask;
                    actID = actID+1;
                end
                
                %precedence
                for l = 1:length(name)-1
                    newPrec = ActivityPrecedence(name(l), name(l+1));
                    newTask = newTask.addPrecedence(newPrec);
                end
                
                %reply-entry: For entry-phase-activities, phase-1 activities reply implicitly
                % Find the last phase-1 activity and set it as the reply activity
                if ~isempty(name) && ~isempty(name{1})
                    % The last phase-1 activity (name{1}) should reply to the entry
                    newEntry.replyActivity{end+1} = name{1};
                end
            end
            
            entries{end+1,1} = newEntry.name;
            entries{end,2} = taskID;
            entries{end,3} = procID;
            newTask = newTask.addEntry(newEntry);
            newEntry.parent = newTask;
            entryID = entryID+1;
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
                    name = char(actElement.getAttribute('name'));
                    hostDemandMean = str2double(char(actElement.getAttribute('host-demand-mean')));
                    hostDemandSCV = str2double(char(actElement.getAttribute('host-demand-cvsq')));
                    if isnan(hostDemandSCV)
                        hostDemandSCV = 1.0;
                    end
                    if hostDemandMean <= 0.0
                        hostDemand = Immediate.getInstance();
                    else
                        if hostDemandSCV <= 0.0
                            hostDemand = Det(hostDemandMean);
                        elseif hostDemandSCV < 1.0
                            hostDemand = APH.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                        elseif hostDemandSCV == 1.0
                            hostDemand = Exp.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                        else
                            hostDemand = HyperExp.fitMeanAndSCV(hostDemandMean, hostDemandSCV);
                        end
                    end
                    boundToEntry = char(actElement.getAttribute('bound-to-entry'));
                    callOrder = char(actElement.getAttribute('call-order'));
                    newAct = Activity(myLN, name, hostDemand, boundToEntry, callOrder);

                    % Parse activity think-time
                    actThinkTimeMean = str2double(char(actElement.getAttribute('think-time')));
                    if ~isnan(actThinkTimeMean) && actThinkTimeMean > 0.0
                        newAct.setThinkTime(actThinkTimeMean);
                    end

                    actObj{end+1,1} = newAct;

                    %synch-call
                    synchCalls = actElement.getElementsByTagName('synch-call');
                    for m = 0:synchCalls.getLength()-1
                        callElement = synchCalls.item(m);
                        dest = char(callElement.getAttribute('dest'));
                        mean = str2double(char(callElement.getAttribute('calls-mean')));
                        newAct = newAct.synchCall(dest,mean);
                    end

                    %asynch-call
                    asynchCalls = actElement.getElementsByTagName('asynch-call');
                    for m = 0:asynchCalls.getLength()-1
                        callElement = asynchCalls.item(m);
                        dest = char(callElement.getAttribute('dest'));
                        mean = str2double(char(callElement.getAttribute('calls-mean')));
                        newAct = newAct.asynchCall(dest,mean);
                    end

                    activities{end+1,1} = newAct.name;
                    activities{end,2} = taskID;
                    activities{end,3} = procID;
                    newTask = newTask.addActivity(newAct);
                    newAct.parent = newTask;
                    actID = actID+1;
                end
            end
            
            %precedence
            precList = taskActsElement.getElementsByTagName('precedence');
            for l = 0:precList.getLength()-1
                precElement = precList.item(l);
                
                %pre
                preTypes = {ActivityPrecedenceType.PRE_SEQ,ActivityPrecedenceType.PRE_AND,ActivityPrecedenceType.PRE_OR};
                for m = 1:length(preTypes)
                    preType = preTypes{m};
                    preList = precElement.getElementsByTagName(ActivityPrecedenceType.toText(preType));
                    if preList.getLength() > 0
                        break
                    end
                end
                preElement = preList.item(0);
                preParams = [];
                preActList = preElement.getElementsByTagName('activity');

                % assumes that preType is a numeric precedence ID
                % if strcmp(preType,ActivityPrecedenceType.toText(ActivityPrecedenceType.PRE_OR))
                if preType == ActivityPrecedenceType.PRE_OR
                    preActs = cell(preActList.getLength(),1);
                    preParams = zeros(postActList.getLength(),1);
                    for m = 0:preActList.getLength()-1
                        preActElement = preActList.item(m);
                        preActs{m+1} = char(preActElement.getAttribute('name'));
                        preParams(m+1) = str2double(char(preActElement.getAttribute('prob')));
                    end
                % elseif strcmp(preType,ActivityPrecedenceType.toText(ActivityPrecedenceType.PRE_AND))
                elseif preType == ActivityPrecedenceType.PRE_AND
                    preActs = cell(preActList.getLength(),1);
                    for m = 0:preActList.getLength()-1
                        preActElement = preActList.item(m);
                        preActs{m+1} = char(preActElement.getAttribute('name'));
                    end
                    preParams = str2double(char(preElement.getAttribute('quorum')));
                else % simple PRE
                    preActs = cell(1,1);
                    preActElement = preActList.item(0);
                    preActs{1} = char(preActElement.getAttribute('name'));
                end
                if isnan(preParams)
                    preParams = [];
                end
                
                %post
                postTypes = {ActivityPrecedenceType.POST_SEQ, ActivityPrecedenceType.POST_AND, ActivityPrecedenceType.POST_OR, ActivityPrecedenceType.POST_LOOP, ActivityPrecedenceType.POST_CACHE};
                for m = 1:length(postTypes)
                    postType = postTypes{m};
                    postList = precElement.getElementsByTagName(ActivityPrecedenceType.toText(postType));
                    if postList.getLength() > 0
                        break
                    end
                end
        
                postElement = postList.item(0);
                postActList = postElement.getElementsByTagName('activity');
                
                % assumes that postType is a numeric precedence ID
                % if strcmp(postType,ActivityPrecedenceType.toText(ActivityPrecedenceType.POST_OR))
                if postType == ActivityPrecedenceType.POST_OR
                    postActs = cell(postActList.getLength(),1);
                    postParams = zeros(postActList.getLength(),1);
                    for m = 0:postActList.getLength()-1
                        postActElement = postActList.item(m);
                        postActs{m+1} = char(postActElement.getAttribute('name'));
                        postParams(m+1) = str2double(char(postActElement.getAttribute('prob')));
                    end
                %elseif strcmp(postType,ActivityPrecedenceType.toText(ActivityPrecedenceType.POST_LOOP))
                elseif postType == ActivityPrecedenceType.POST_LOOP
                    postActs = cell(postActList.getLength()+1,1);
                    postParams = zeros(postActList.getLength(),1);
                    for m = 0:postActList.getLength()-1
                        postActElement = postActList.item(m);
                        postActs{m+1} = char(postActElement.getAttribute('name'));
                        postParams(m+1) = str2double(char(postActElement.getAttribute('count')));
         
                    end
                    postActs{end} = char(postElement.getAttribute('end'));                    
                else
                    postActs = cell(postActList.getLength(),1);
                    postParams = [];
                    for m = 0:postActList.getLength()-1
                        postActElement = postActList.item(m);
                        postActs{m+1} = char(postActElement.getAttribute('name'));
                    end
                end
                newPrec = ActivityPrecedence(preActs, postActs, preType, postType, preParams, postParams);
                newTask = newTask.addPrecedence(newPrec);
            end
            
            %reply-entry
            replyList = taskActsElement.getElementsByTagName('reply-entry');
            for l = 0:replyList.getLength()-1
                replyElement = replyList.item(l);
                replyName = char(replyElement.getAttribute('name'));
                replyIdx = findstring(entries(:,1), replyName);
                replyActList = replyElement.getElementsByTagName('reply-activity');
                for m = 0:replyActList.getLength()-1
                    replyActElement = replyActList.item(m);
                    replyActName = char(replyActElement.getAttribute('name'));
                    entryObj{replyIdx}.replyActivity{end+1} = replyActName;
                end
            end
        end
        
        tasks{end+1,1} = newTask.name;
        tasks{end,2} = procID;
        newProc = newProc.addTask(newTask);
        taskID = taskID+1;
    end
    
    hosts{end+1,1} = newProc.name;
    procID = procID+1;
end
end
