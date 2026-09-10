function myLN = parseXML(filename, verbose)
% MYLN = PARSEXML(FILENAME, VERBOSE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


import javax.xml.parsers.DocumentBuilderFactory; %#ok<JAPIEXT629>
import javax.xml.parsers.DocumentBuilder;        %#ok<JAPIEXT629>
import org.w3c.dom.Document;                     %#ok<JAPIEXT161>
import org.w3c.dom.NodeList;                     %#ok<JAPIEXT161>
import org.w3c.dom.Node;                         %#ok<JAPIEXT161>
import org.w3c.dom.Element;                      %#ok<JAPIEXT161>
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
validateInputModel(doc);
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
        % LINE dialect (.lqnx cache/setup extension): a <cache> child makes the
        % task a CacheTask, a <setup> or <delay-off> child makes it a SetupTask.
        % Both constructors take neither a think time nor a reply, so the think
        % time is applied after construction rather than through the Task one.
        cacheList = child_elements(taskElement, 'cache');
        setupList = child_elements(taskElement, 'setup');
        delayOffList = child_elements(taskElement, 'delay-off');
        if ~isempty(cacheList)
            cacheElement = cacheList{1};
            nitems = str2double(char(cacheElement.getAttribute('items')));
            levelList = child_elements(cacheElement, 'level');
            caps = zeros(1, numel(levelList));
            for lv = 1:numel(levelList)
                caps(lv) = str2double(char(levelList{lv}.getAttribute('capacity')));
            end
            if isempty(levelList)
                line_error(mfilename, sprintf(['Task "%s" declares a <cache> with no <level> child: ' ...
                    'the cache capacity is one <level capacity="..."/> per cache list.'], name));
            end
            replStr = char(cacheElement.getAttribute('replacement'));
            newTask = CacheTask(myLN, name, nitems, caps, repl_from_lqnx(replStr, name), ...
                multiplicity, SchedStrategy.fromText(scheduling));
            retrievalStr = char(cacheElement.getAttribute('retrieval'));
            if strcmpi(retrievalStr, 'true') || strcmp(retrievalStr, '1')
                newTask.setRetrieval(true);
            end
            newTask.setThinkTime(thinkTime);
        elseif ~isempty(setupList) || ~isempty(delayOffList)
            newTask = SetupTask(myLN, name, multiplicity, SchedStrategy.fromText(scheduling));
            newTask.setThinkTime(thinkTime);
        else
            newTask = Task(myLN, name, multiplicity, SchedStrategy.fromText(scheduling), thinkTime);
        end
        if ~isempty(setupList)
            newTask.setSetupTime(time_from_element(setupList{1}));
        end
        if ~isempty(delayOffList)
            newTask.setDelayOffTime(time_from_element(delayOffList{1}));
        end
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

        % Parse fan-out elements if present (used for replication load distribution)
        fanOutList = taskElement.getElementsByTagName('fan-out');
        for fo = 0:fanOutList.getLength()-1
            fanOutElement = fanOutList.item(fo);
            dest = char(fanOutElement.getAttribute('dest'));
            valueStr = char(fanOutElement.getAttribute('value'));
            if ~isempty(dest) && ~isempty(valueStr)
                value = str2double(valueStr);
                if ~isnan(value)
                    newTask.setFanOut(dest, value);
                end
            end
        end

        taskObj{end+1,1} = newTask;

        entryList = taskElement.getElementsByTagName('entry');
        for k = 0:entryList.getLength()-1
            %Element - Entry
            entryElement = entryList.item(k);
            name = char(entryElement.getAttribute('name'));
            % LINE dialect: an <item-entry> child makes the entry an ItemEntry.
            itemList = child_elements(entryElement, 'item-entry');
            if isempty(itemList)
                newEntry = Entry(myLN, name);
            else
                itemElement = itemList{1};
                cardinality = str2double(char(itemElement.getAttribute('cardinality')));
                popularity = popularity_from_element(itemElement, cardinality, name);
                newEntry = ItemEntry(myLN, name, cardinality, popularity);
            end
            openArrivalRate = str2double(char(entryElement.getAttribute('open-arrival-rate')));
            if ~isnan(openArrivalRate) && openArrivalRate > 0
                newEntry.setArrival(Exp.fitMean(1/openArrivalRate));
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
                    newAct.setPhase(phase);  % validated, not assigned past the guard

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

                    %call-group (LINE dialect)
                    newAct = parse_call_groups(actElement, newAct);

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

                    %call-group (LINE dialect)
                    newAct = parse_call_groups(actElement, newAct);

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
                % The post side is minOccurs="0" in lqn-core.xsd: a precedence
                % that carries only a pre element declares a TERMINAL activity
                % and no successor, so it contributes no edge and is skipped.
                postTypes = {ActivityPrecedenceType.POST_SEQ, ActivityPrecedenceType.POST_AND, ActivityPrecedenceType.POST_OR, ActivityPrecedenceType.POST_LOOP, ActivityPrecedenceType.POST_CACHE};
                hasPost = false;
                for m = 1:length(postTypes)
                    postType = postTypes{m};
                    postList = precElement.getElementsByTagName(ActivityPrecedenceType.toText(postType));
                    if postList.getLength() > 0
                        hasPost = true;
                        break
                    end
                end
                if ~hasPost
                    continue
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
                elseif postType == ActivityPrecedenceType.POST_CACHE
                    % LINE dialect: the branch is NAMED by cache-result, so a
                    % reordered list cannot swap hit for miss. A file written
                    % before the attribute existed carries none, and document
                    % order (hit first) is then the answer, which is what it
                    % meant when it was written.
                    postActs = cell(postActList.getLength(),1);
                    postParams = [];
                    results = cell(postActList.getLength(),1);
                    for m = 0:postActList.getLength()-1
                        postActElement = postActList.item(m);
                        postActs{m+1} = char(postActElement.getAttribute('name'));
                        results{m+1} = lower(char(postActElement.getAttribute('cache-result')));
                    end
                    named = ~cellfun(@isempty, results);
                    if any(named)
                        hitIdx = find(strcmp(results, 'hit'), 1);
                        missIdx = find(strcmp(results, 'miss'), 1);
                        if ~isempty(hitIdx) && ~isempty(missIdx)
                            rest = setdiff(1:numel(postActs), [hitIdx, missIdx], 'stable');
                            postActs = postActs([hitIdx, missIdx, rest]);
                        end
                    end
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

function validateInputModel(doc)
% VALIDATEINPUTMODEL Reject a structurally inconsistent LQN document
%
% Run on the parsed document before any object is built, so that a defective
% input is named at its source instead of surfacing as a downstream failure.
% The same checks, in the same order and with the same messages, are applied
% by the JAR, Python and C++ readers.

tol = 1e-6;
procNames = cell(0,1);
taskNames = cell(0,1);
entryNames = cell(0,1);
entryOwner = cell(0,1); % task owning entryNames{k}
isRefEntry = false(0,1);
callDests = cell(0,1);
replyEntries = cell(0,1);
hasRefTask = false;
hasOpenArrival = false;

procList = doc.getElementsByTagName('processor');
for i = 0:procList.getLength()-1
    procElement = procList.item(i);
    procName = char(procElement.getAttribute('name'));
    if any(strcmp(procNames, procName))
        line_error(mfilename, sprintf('Duplicate processor name "%s".', procName));
    end
    procNames{end+1,1} = procName; %#ok<AGROW>

    taskList = procElement.getElementsByTagName('task');
    for j = 0:taskList.getLength()-1
        taskElement = taskList.item(j);
        taskName = char(taskElement.getAttribute('name'));
        if any(strcmp(taskNames, taskName))
            line_error(mfilename, sprintf('Duplicate task name "%s".', taskName));
        end
        taskNames{end+1,1} = taskName; %#ok<AGROW>
        isRef = strcmpi(char(taskElement.getAttribute('scheduling')), 'ref');
        hasRefTask = hasRefTask || isRef;

        entryList = taskElement.getElementsByTagName('entry');
        if entryList.getLength() == 0
            line_error(mfilename, sprintf('Task "%s" has no entries.', taskName));
        end
        for k = 0:entryList.getLength()-1
            entryElement = entryList.item(k);
            entryName = char(entryElement.getAttribute('name'));
            if any(strcmp(entryNames, entryName))
                line_error(mfilename, sprintf('Duplicate entry name "%s".', entryName));
            end
            entryNames{end+1,1} = entryName; %#ok<AGROW>
            entryOwner{end+1,1} = taskName; %#ok<AGROW>
            isRefEntry(end+1,1) = isRef; %#ok<AGROW>

            openArrivalRate = str2double(char(entryElement.getAttribute('open-arrival-rate')));
            if ~isnan(openArrivalRate) && openArrivalRate > 0
                hasOpenArrival = true;
                if isRef
                    line_error(mfilename, sprintf('Entry "%s" belongs to reference task "%s" and cannot have open arrivals.', entryName, taskName));
                end
            end

            fwdList = entryElement.getElementsByTagName('forwarding');
            if isRef && fwdList.getLength() > 0
                line_error(mfilename, sprintf('Entry "%s" belongs to reference task "%s" and cannot forward requests.', entryName, taskName));
            end
            fwdTotal = 0.0;
            for fw = 0:fwdList.getLength()-1
                fwdElement = fwdList.item(fw);
                probStr = char(fwdElement.getAttribute('prob'));
                if isempty(probStr)
                    prob = 1.0;
                else
                    prob = str2double(probStr);
                end
                if isnan(prob) || prob < 0.0 || prob > 1.0
                    line_error(mfilename, sprintf('Forwarding from entry "%s" to entry "%s" has an invalid probability of %g.', entryName, char(fwdElement.getAttribute('dest')), prob));
                end
                fwdTotal = fwdTotal + prob;
            end
            if fwdTotal > 1.0 + tol
                line_error(mfilename, sprintf('Entry "%s" has a total forwarding probability of %g.', entryName, fwdTotal));
            end
        end

        % activity names are unique within their task; a name under a pre or post list is a reference, not a declaration
        actNames = cell(0,1);
        actList = taskElement.getElementsByTagName('activity');
        for l = 0:actList.getLength()-1
            actElement = actList.item(l);
            parentTag = char(actElement.getParentNode().getNodeName());
            if ~strcmp(parentTag,'task-activities') && ~strcmp(parentTag,'entry-phase-activities')
                continue
            end
            actName = char(actElement.getAttribute('name'));
            if any(strcmp(actNames, actName))
                line_error(mfilename, sprintf('Duplicate activity name "%s" in task "%s".', actName, taskName));
            end
            actNames{end+1,1} = actName; %#ok<AGROW>
        end

        callList = taskElement.getElementsByTagName('synch-call');
        for m = 0:callList.getLength()-1
            callDests{end+1,1} = char(callList.item(m).getAttribute('dest')); %#ok<AGROW>
        end
        callList = taskElement.getElementsByTagName('asynch-call');
        for m = 0:callList.getLength()-1
            callDests{end+1,1} = char(callList.item(m).getAttribute('dest')); %#ok<AGROW>
        end
        fwdList = taskElement.getElementsByTagName('forwarding');
        for fw = 0:fwdList.getLength()-1
            callDests{end+1,1} = char(fwdList.item(fw).getAttribute('dest')); %#ok<AGROW>
        end

        orList = taskElement.getElementsByTagName('post-OR');
        for l = 0:orList.getLength()-1
            branchList = orList.item(l).getElementsByTagName('activity');
            branchTotal = 0.0;
            for m = 0:branchList.getLength()-1
                branchElement = branchList.item(m);
                probStr = char(branchElement.getAttribute('prob'));
                if isempty(probStr)
                    prob = 1.0;
                else
                    prob = str2double(probStr);
                end
                if isnan(prob) || prob < 0.0 || prob > 1.0
                    line_error(mfilename, sprintf('Activity "%s" in task "%s" has an invalid branch probability of %g.', char(branchElement.getAttribute('name')), taskName, prob));
                end
                branchTotal = branchTotal + prob;
            end
            if abs(branchTotal - 1.0) > tol
                line_error(mfilename, sprintf('Branch probabilities of an OR-fork in task "%s" sum to %g instead of 1.', taskName, branchTotal));
            end
        end

        replyList = taskElement.getElementsByTagName('reply-entry');
        for l = 0:replyList.getLength()-1
            replyEntries{end+1,1} = char(replyList.item(l).getAttribute('name')); %#ok<AGROW>
        end
    end
end

for c = 1:length(callDests)
    idx = find(strcmp(entryNames, callDests{c}), 1);
    if ~isempty(idx) && isRefEntry(idx)
        line_error(mfilename, sprintf('Entry "%s" belongs to reference task "%s" and cannot receive requests.', entryNames{idx}, entryOwner{idx}));
    end
end

for r = 1:length(replyEntries)
    idx = find(strcmp(entryNames, replyEntries{r}), 1);
    if ~isempty(idx) && isRefEntry(idx)
        line_error(mfilename, sprintf('Entry "%s" belongs to reference task "%s" and cannot be replied to.', entryNames{idx}, entryOwner{idx}));
    end
end

if ~hasRefTask && ~hasOpenArrival
    line_error(mfilename, 'The model has no reference task and no open arrivals.');
end
end


function kids = child_elements(parent, tagName)
% KIDS = CHILD_ELEMENTS(PARENT, TAGNAME)
% The DIRECT children of PARENT named TAGNAME, as a cell array.
%
% getElementsByTagName searches every DESCENDANT, which is wrong for an element
% that may legally nest: a <cache> found under a task would otherwise be claimed
% by an enclosing element as well. Walking the child list keeps each element with
% the node that declares it, and an unknown child is simply not matched, so a
% reader meeting a future element does not fail.
kids = {};
nodes = parent.getChildNodes();
for i = 0:nodes.getLength()-1
    node = nodes.item(i);
    if node.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE && ...
            strcmp(char(node.getNodeName()), tagName)
        kids{end+1} = node; %#ok<AGROW>
    end
end
end


function act = parse_call_groups(actElement, act)
% ACT = PARSE_CALL_GROUPS(ACTELEMENT, ACT)
% LINE dialect <call-group>: which of the synch-calls parsed above one
% dispatcher issues, and under which strategy. The member calls are ordinary
% synch-call elements and have already been read, so only the grouping is
% recorded; issuing them again here would double the call rate.
grps = child_elements(actElement, 'call-group');
for g = 1:numel(grps)
    strategy = callgroup_from_lqnx(char(grps{g}.getAttribute('strategy')), act.name);
    destElems = child_elements(grps{g}, 'dest');
    dests = cell(1, numel(destElems));
    for d = 1:numel(destElems)
        dests{d} = char(destElems{d}.getAttribute('name'));
    end
    act = act.recordCallGroup(strategy, dests);
end
end


function strategy = callgroup_from_lqnx(name, actName)
% STRATEGY = CALLGROUP_FROM_LQNX(NAME, ACTNAME)
% Wire enum name -> RoutingStrategy id, the inverse of writeXML's
% callGroupStrategyName.
switch upper(strtrim(name))
    case 'RROBIN', strategy = RoutingStrategy.RROBIN;
    case 'JSQ',    strategy = RoutingStrategy.JSQ;
    otherwise
        line_error(mfilename, sprintf(['Activity "%s" declares a call group with an ' ...
            'unrecognized strategy "%s"; the dialect spells them RROBIN and JSQ.'], ...
            actName, name));
end
end


function id = repl_from_lqnx(name, taskName)
% ID = REPL_FROM_LQNX(NAME, TASKNAME)
% Wire enum name -> ReplacementStrategy id, the inverse of writeXML's
% repl_to_lqnx and of linemodel_save's repl_to_str.
switch upper(strtrim(name))
    case 'LRU',   id = ReplacementStrategy.LRU;
    case 'FIFO',  id = ReplacementStrategy.FIFO;
    case 'RR',    id = ReplacementStrategy.RR;
    case 'SFIFO', id = ReplacementStrategy.SFIFO;
    case 'HLRU',  id = ReplacementStrategy.HLRU;
    case 'CLIMB', id = ReplacementStrategy.CLIMB;
    case 'QLRU',  id = ReplacementStrategy.QLRU;
    otherwise
        line_error(mfilename, sprintf(['Task "%s" declares an unrecognized cache replacement ' ...
            'strategy "%s"; the dialect spells them RR, FIFO, LRU, SFIFO, HLRU, CLIMB, QLRU.'], ...
            taskName, name));
end
end


function dist = time_from_element(element)
% DIST = TIME_FROM_ELEMENT(ELEMENT)
% A <setup>/<delay-off> mean and SCV back into a distribution, taking the same
% family the setter itself would: an SCV of one is the Exp the numeric
% setSetupTime builds, anything else needs a two-moment fit.
mean = str2double(char(element.getAttribute('mean')));
scv = str2double(char(element.getAttribute('scv')));
if isnan(scv)
    scv = 1.0;
end
if isnan(mean) || mean <= GlobalConstants.FineTol
    dist = Immediate.getInstance();
elseif abs(scv - 1.0) <= GlobalConstants.FineTol
    dist = Exp.fitMean(mean);
else
    dist = APH.fitMeanAndSCV(mean, scv);
end
end


function dist = popularity_from_element(itemElement, cardinality, entryName)
% DIST = POPULARITY_FROM_ELEMENT(ITEMELEMENT, CARDINALITY, ENTRYNAME)
% The inverse of writeXML's popularity_element. The parameter list is in
% CONSTRUCTOR ORDER; a DiscreteSampler's is split on the cardinality already
% declared on the parent <item-entry>, so n values are p over the default
% support 1..n and 2n are p followed by x.
popList = child_elements(itemElement, 'access-popularity');
if isempty(popList)
    % No popularity was written, which the dialect allows. Uniform access over
    % the declared cardinality is the only reading that is not a guess.
    dist = DiscreteSampler(ones(1, cardinality) / cardinality);
    return
end
popElement = popList{1};
cls = char(popElement.getAttribute('name'));
paramList = child_elements(popElement, 'parameter');
vals = zeros(1, numel(paramList));
for i = 1:numel(paramList)
    vals(i) = str2double(char(paramList{i}.getAttribute('value')));
end
switch cls
    case 'DiscreteSampler'
        n = cardinality;
        if numel(vals) == 2 * n
            dist = DiscreteSampler(vals(1:n), vals(n+1:end));
        elseif numel(vals) == n
            dist = DiscreteSampler(vals);
        else
            line_error(mfilename, sprintf(['Entry "%s" declares a DiscreteSampler popularity with ' ...
                '%d parameters, which is neither the cardinality %d (p alone) nor twice it ' ...
                '(p followed by the support x).'], entryName, numel(vals), n));
        end
    case 'Zipf'
        if numel(vals) ~= 2
            line_error(mfilename, sprintf(['Entry "%s" declares a Zipf popularity with %d ' ...
                'parameters; Zipf(s, n) takes two.'], entryName, numel(vals)));
        end
        dist = Zipf(vals(1), vals(2));
    otherwise
        line_error(mfilename, sprintf(['Entry "%s" declares a popularity distribution of class ' ...
            '''%s'', which this reader cannot rebuild from a flat parameter list. Extend ' ...
            'popularity_element/popularity_from_element in writeXML.m and parseXML.m together.'], ...
            entryName, cls));
end
end
