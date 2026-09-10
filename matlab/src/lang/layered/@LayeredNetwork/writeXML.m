function writeXML(self,filename,useAbstractNames)
% WRITEXML(SELF,FILENAME,USEAN)
%
% USEAN: if true replaces node names with abstract names eg E1, T1, ...
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<3
    useAbstractNames=false;
end

% The LQN schema has no attribute for a service-rate dependence or for a
% compatibility graph over a server's operands, so neither can survive the file;
% warn rather than let a reader see a plain homogeneous server. A pool
% declaration is named separately because it is a structural statement about the
% servers, not a scaling on top of them, and a reader who loses it silently gets
% a DIFFERENT system of the same total multiplicity.
for p = 1:length(self.hosts)
    if self.hosts{p}.hasServerPools()
        line_warning(mfilename,'Host %s declares heterogeneous server pools with a class-compatibility graph, which the LQN XML schema cannot represent: they are omitted from %s and the host reads back as a homogeneous multiserver.', self.hosts{p}.name, filename);
    elseif self.hosts{p}.hasRateDependence()
        line_warning(mfilename,'Host %s declares a service-rate dependence, which the LQN XML schema cannot represent: it is omitted from %s.', self.hosts{p}.name, filename);
    end
end
for t = 1:length(self.tasks)
    if self.tasks{t}.hasServerPools()
        line_warning(mfilename,'Task %s declares heterogeneous server pools with a class-compatibility graph, which the LQN XML schema cannot represent: they are omitted from %s and the task reads back as a homogeneous multiserver.', self.tasks{t}.name, filename);
    elseif self.tasks{t}.hasRateDependence()
        line_warning(mfilename,'Task %s declares a service-rate dependence, which the LQN XML schema cannot represent: it is omitted from %s.', self.tasks{t}.name, filename);
    end
end

nodeHashMap = configureDictionary('string','cell');

tctr = 0;
ectr = 0;
actr = 0;
if useAbstractNames
    for p = 1:length(self.hosts)
        curProc = self.hosts{p};
        nodeHashMap{curProc.name}=sprintf('P%d',p);
        for t=1:length(curProc.tasks)
            curTask = curProc.tasks(t);
            tctr = tctr + 1;
            nodeHashMap{curTask.name}=sprintf('T%d',tctr);
            for e=1:length(curTask.entries)
                curEntry = curTask.entries(e);
                ectr = ectr + 1;
                nodeHashMap{curEntry.name}=sprintf('E%d',ectr);
            end
            for a=1:length(curTask.activities)
                curAct = curTask.activities(a);
                actr = actr + 1;
                nodeHashMap{curAct.name} = sprintf('A%d',actr);
            end
        end
    end
else
    for p = 1:length(self.hosts)
        curProc = self.hosts{p};
        nodeHashMap{curProc.name}=curProc.name;
        for t=1:length(curProc.tasks)
            curTask = curProc.tasks(t);
            nodeHashMap{curTask.name}=curTask.name;
            for e=1:length(curTask.entries)
                curEntry = curTask.entries(e);
                nodeHashMap{curEntry.name}=curEntry.name;
            end
            for a=1:length(curTask.activities)
                curAct = curTask.activities(a);
                nodeHashMap{curAct.name}=curAct.name;
            end
        end
    end
end

import javax.xml.parsers.DocumentBuilderFactory; %#ok<JAPIEXT629>
import javax.xml.parsers.DocumentBuilder;        %#ok<JAPIEXT629>
import org.w3c.dom.Document;                     %#ok<JAPIEXT161>
import org.w3c.dom.NodeList;                     %#ok<JAPIEXT161>
import org.w3c.dom.Node;                         %#ok<JAPIEXT161>
import org.w3c.dom.Element;                      %#ok<JAPIEXT161>
import java.io.File;
import javax.xml.transform.Transformer;          %#ok<JAPIEXT216>
import javax.xml.transform.TransformerFactory;   %#ok<JAPIEXT216>
import javax.xml.transform.dom.DOMSource;        %#ok<JAPIEXT216>
import javax.xml.transform.stream.StreamResult;  %#ok<JAPIEXT216>
import javax.xml.transform.OutputKeys;           %#ok<JAPIEXT216>

precision = '%10.15e'; %precision for doubles
docFactory = DocumentBuilderFactory.newInstance();
docBuilder = docFactory.newDocumentBuilder();
doc = docBuilder.newDocument();

%Root Element
rootElement = doc.createElement('lqn-model');
doc.appendChild(rootElement);
rootElement.setAttribute('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance');
rootElement.setAttribute('xsi:noNamespaceSchemaLocation', 'lqn.xsd');
rootElement.setAttribute('name', getName(self));

for p = 1:length(self.hosts)
    %processor
    curProc = self.hosts{p};
    procElement = doc.createElement('processor');
    rootElement.appendChild(procElement);
    procElement.setAttribute('name', nodeHashMap{curProc.name});
    procElement.setAttribute('scheduling', lqnx_sched(curProc.scheduling));
    if curProc.replication>1
        procElement.setAttribute('replication', num2str(curProc.replication));
    end
    if SchedStrategy.fromText(curProc.scheduling) ~= SchedStrategy.INF
        mult = num2str(curProc.multiplicity);
        if isinf(mult), mult=1; end
        procElement.setAttribute('multiplicity', mult);
    end
    if SchedStrategy.fromText(curProc.scheduling) == SchedStrategy.PS || SchedStrategy.fromText(curProc.scheduling) == SchedStrategy.PSPRIO
        procElement.setAttribute('quantum', lqnx_num(curProc.quantum));
    end
    procElement.setAttribute('speed-factor', lqnx_num(curProc.speedFactor));
    for t=1:length(curProc.tasks)
        curTask = curProc.tasks(t);
        taskElement = doc.createElement('task');
        procElement.appendChild(taskElement);
        taskElement.setAttribute('name', nodeHashMap{curTask.name});
        taskElement.setAttribute('scheduling', lqnx_sched(curTask.scheduling));
        % parseXML reads this back into Task.priority, so omitting it would drop
        % the priority of a prioritised task on a write/read round trip
        if curTask.priority ~= 0
            taskElement.setAttribute('priority', num2str(curTask.priority));
        end
        if curTask.replication>1
            taskElement.setAttribute('replication',  num2str(curTask.replication));
        end
        if  SchedStrategy.fromText(curTask.scheduling) ~= SchedStrategy.INF
            taskElement.setAttribute('multiplicity', num2str(curTask.multiplicity));
        end
        % think-time may only be written for a reference task: lqns rejects the
        % file outright ('Task "X" is not a reference task; it cannot have think
        % time'). LINE does give a non-reference task's think time to its
        % callers, so that value cannot survive this format and the loss is
        % reported rather than left silent. Bridges that need it (PYLINE
        % lang='python') reapply it after parseXML.
        if SchedStrategy.fromText(curTask.scheduling) == SchedStrategy.REF
            taskElement.setAttribute('think-time', lqnx_num(curTask.thinkTimeMean));
        elseif ~isempty(curTask.thinkTimeMean) && curTask.thinkTimeMean > 0
            line_warning(mfilename, sprintf(['Task %s is not a reference task, so its think time (%g) ', ...
                'is not written: the LQN XML schema accepts think-time on reference tasks only.'], ...
                curTask.name, curTask.thinkTimeMean));
        end
        % LINE dialect (see the .lqnx cache/setup extension): a CacheTask, an
        % ItemEntry and a task's setup/delay-off times have no element in
        % lqn-core.xsd, and writing the file without them emitted a VALID LQN
        % describing a DIFFERENT model -- a plain FCFS task with no cache and no
        % setup, which every reader then solved correctly for the wrong thing.
        % These elements carry the same fields the JSON interchange already
        % carries (linemodel_save: totalItems, cacheCapacity, replacementStrategy,
        % accessProb, setupTime, delayOffTime), so the two transports agree.
        if isa(curTask, 'CacheTask')
            cacheElement = doc.createElement('cache');
            taskElement.appendChild(cacheElement);
            cacheElement.setAttribute('items', num2str(curTask.items));
            cacheElement.setAttribute('replacement', repl_to_lqnx(curTask.replacestrategy));
            if curTask.hasRetrieval()
                cacheElement.setAttribute('retrieval', 'true');
            end
            % ONE <level> PER CACHE LIST, IN ORDER: itemLevelCap is an array and a
            % multi-list cache is the normal case, so a scalar attribute would
            % silently merge the lists into one.
            caps = curTask.itemLevelCap;
            for lv = 1:numel(caps)
                levelElement = doc.createElement('level');
                cacheElement.appendChild(levelElement);
                levelElement.setAttribute('capacity', num2str(caps(lv)));
            end
        end
        if isprop(curTask, 'setupTimeMean') && ~isempty(curTask.setupTimeMean) ...
                && curTask.setupTimeMean > GlobalConstants.FineTol
            setupElement = doc.createElement('setup');
            taskElement.appendChild(setupElement);
            setupElement.setAttribute('mean', lqnx_num(curTask.setupTimeMean));
            setupElement.setAttribute('scv', lqnx_num(curTask.setupTimeSCV));
        end
        if isprop(curTask, 'delayOffTimeMean') && ~isempty(curTask.delayOffTimeMean) ...
                && curTask.delayOffTimeMean > GlobalConstants.FineTol
            delayOffElement = doc.createElement('delay-off');
            taskElement.appendChild(delayOffElement);
            delayOffElement.setAttribute('mean', lqnx_num(curTask.delayOffTimeMean));
            delayOffElement.setAttribute('scv', lqnx_num(curTask.delayOffTimeSCV));
        end
        % fan-out/fan-in are parsed by parseXML, by the JAR and by the cpp
        % lqn_reader, and were written by no codebase, so a replicated model
        % lost its call multiplicities on every round trip. lqn-core.xsd
        % (TaskType) places them before the entries.
        for fo=1:length(curTask.fanOutDest)
            fanOutElement = doc.createElement('fan-out');
            taskElement.appendChild(fanOutElement);
            fanOutElement.setAttribute('dest', nodeHashMap{curTask.fanOutDest{fo}});
            fanOutElement.setAttribute('value', lqnx_num(curTask.fanOutValue(fo)));
        end
        if ~isempty(curTask.fanInSource) && curTask.fanInValue > 0
            fanInElement = doc.createElement('fan-in');
            taskElement.appendChild(fanInElement);
            fanInElement.setAttribute('source', nodeHashMap{curTask.fanInSource});
            fanInElement.setAttribute('value', lqnx_num(curTask.fanInValue));
        end
        for e=1:length(curTask.entries)
            curEntry = curTask.entries(e);
            entryElement = doc.createElement('entry');
            taskElement.appendChild(entryElement);
            entryElement.setAttribute('name', nodeHashMap{curEntry.name});
            entryElement.setAttribute('type', 'NONE');
            if ~isempty(curEntry.arrival) && isa(curEntry.arrival, 'Distribution') && ...
                    isfinite(curEntry.arrival.getMean()) && curEntry.arrival.getMean() > GlobalConstants.FineTol
                entryElement.setAttribute('open-arrival-rate', lqnx_num(1 / curEntry.arrival.getMean()));
            end
            % LINE dialect: the item reference stream of an ItemEntry. Its
            % presence is what makes the entry an ItemEntry on read.
            if isa(curEntry, 'ItemEntry')
                itemElement = doc.createElement('item-entry');
                entryElement.appendChild(itemElement);
                itemElement.setAttribute('cardinality', num2str(curEntry.cardinality));
                if ~isempty(curEntry.popularity) && isa(curEntry.popularity, 'Distribution')
                    popElement = popularity_element(doc, curEntry.popularity, ...
                        curEntry.cardinality, curEntry.name);
                    itemElement.appendChild(popElement);
                end
            end

            % Write forwarding calls
            for fw=1:length(curEntry.forwardingDests)
                fwdElement = doc.createElement('forwarding');
                entryElement.appendChild(fwdElement);
                fwdElement.setAttribute('dest', nodeHashMap{curEntry.forwardingDests{fw}});
                fwdElement.setAttribute('prob', lqnx_num(curEntry.forwardingProbs(fw)));
            end
        end
        taskActsElement = doc.createElement('task-activities');
        taskElement.appendChild(taskActsElement);
        for a=1:length(curTask.activities)
            curAct = curTask.activities(a);
            actElement = doc.createElement('activity');
            taskActsElement.appendChild(actElement);
            actElement.setAttribute('host-demand-mean', lqnx_num(curAct.hostDemandMean));
            actElement.setAttribute('host-demand-cvsq', lqnx_num(curAct.hostDemandSCV));
            if ~isempty(curAct.boundToEntry)
                actElement.setAttribute('bound-to-entry', nodeHashMap{curAct.boundToEntry});
            end
            actElement.setAttribute('call-order', curAct.callOrder);
            actElement.setAttribute('name', nodeHashMap{curAct.name});
            if curAct.thinkTimeMean > GlobalConstants.FineTol
                actElement.setAttribute('think-time', lqnx_num(curAct.thinkTimeMean));
            end

            for sc=1:length(curAct.syncCallDests)
                syncCallElement = doc.createElement('synch-call');
                actElement.appendChild(syncCallElement);
                syncCallElement.setAttribute('dest', nodeHashMap{curAct.syncCallDests{sc}});
                syncCallElement.setAttribute('calls-mean', lqnx_num(curAct.syncCallMeans(sc)));
            end
            for ac=1:length(curAct.asyncCallDests)
                asyncCallElement = doc.createElement('asynch-call');
                actElement.appendChild(asyncCallElement);
                asyncCallElement.setAttribute('dest', nodeHashMap{curAct.asyncCallDests{ac}});
                asyncCallElement.setAttribute('calls-mean', lqnx_num(curAct.asyncCallMeans(ac)));
            end
            % LINE dialect: a routed call group names which of the synch-calls
            % above one dispatcher issues, and under which strategy. The member
            % calls stay ordinary synch-calls, so a reader that ignores this
            % element still gets the same aggregate call means -- which is
            % exactly what lqns and lqsim, having no dispatcher, should see.
            for g=1:length(curAct.syncCallGroups)
                grp = curAct.syncCallGroups{g};
                grpElement = doc.createElement('call-group');
                actElement.appendChild(grpElement);
                grpElement.setAttribute('strategy', callGroupStrategyName(grp.strategy));
                for d=1:length(grp.dests)
                    destElement = doc.createElement('dest');
                    grpElement.appendChild(destElement);
                    destElement.setAttribute('name', nodeHashMap{grp.dests{d}});
                end
            end
        end
        for ap=1:length(curTask.precedences)
            curActPrec = curTask.precedences(ap);
            actPrecElement = doc.createElement('precedence');
            taskActsElement.appendChild(actPrecElement);

            preElement = doc.createElement(ActivityPrecedenceType.toText(curActPrec.preType));
            actPrecElement.appendChild(preElement);
            if curActPrec.preType== ActivityPrecedenceType.PRE_AND && ~isempty(curActPrec.preParams)
                preElement.setAttribute('quorum', num2str(curActPrec.preParams(1)));
            end
            for pra = 1:length(curActPrec.preActs)
                preActElement = doc.createElement('activity');
                preElement.appendChild(preActElement);
                preActElement.setAttribute('name', nodeHashMap{curActPrec.preActs{pra}});
            end

            postElement = doc.createElement(ActivityPrecedenceType.toText(curActPrec.postType));
            actPrecElement.appendChild(postElement);
            if curActPrec.postType==ActivityPrecedenceType.POST_OR
                for poa = 1:length(curActPrec.postActs)
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', nodeHashMap{curActPrec.postActs{poa}});
                    postActElement.setAttribute('prob', lqnx_num(curActPrec.postParams(poa)));
                end
            elseif curActPrec.postType==ActivityPrecedenceType.POST_LOOP
                % ActivityPrecedence.Loop carries ONE count for the whole body,
                % which the schema states per looped activity: a body of two or
                % more activities used to index past the end of postParams and
                % the model could not be written out at all
                for poa = 1:length(curActPrec.postActs)-1
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', nodeHashMap{curActPrec.postActs{poa}});
                    if isscalar(curActPrec.postParams)
                        loopCount = curActPrec.postParams;
                    else
                        loopCount = curActPrec.postParams(poa);
                    end
                    postActElement.setAttribute('count', num2str(loopCount));
                end
                % the end activity is named like every other one, i.e. through
                % the hash map, which abbreviated mode renames by
                postElement.setAttribute('end', nodeHashMap{curActPrec.postActs{end}});
            elseif curActPrec.postType==ActivityPrecedenceType.POST_CACHE
                % <post-CACHE> was already emitted (ActivityPrecedenceType.toText
                % names it), but its branches were distinguished by POSITION
                % alone. Name the branch, so a reader that reorders the list, or
                % a writer that does, cannot swap hit for miss.
                names = {'hit', 'miss'};
                for poa = 1:length(curActPrec.postActs)
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', nodeHashMap{curActPrec.postActs{poa}});
                    if poa <= numel(names)
                        postActElement.setAttribute('cache-result', names{poa});
                    end
                end
            else
                for poa = 1:length(curActPrec.postActs)
                    postActElement = doc.createElement('activity');
                    postElement.appendChild(postActElement);
                    postActElement.setAttribute('name', nodeHashMap{curActPrec.postActs{poa}});
                end
            end
        end
        if SchedStrategy.fromText(curTask.scheduling) ~= SchedStrategy.REF
            % Get the model structure to check for sync vs async calls
            lsn = self.getStruct();
            for e=1:length(curTask.entries)
                curEntry = curTask.entries(e);

                % Find this entry's index in the model
                eidx = find(cellfun(@(x) strcmp(x.name, curEntry.name), self.entries));
                globalEidx = lsn.eshift + eidx;

                % Check if this entry is called synchronously or is a forwarding target
                % If only called asynchronously (and not a forwarding target), skip generating reply-entry
                isCalledSync = any(lsn.issynccaller(:, globalEidx));

                % Check if this entry is a forwarding target (receives forwarded requests)
                isForwardingTarget = false;
                for cidx = 1:lsn.ncalls
                    if full(lsn.calltype(cidx)) == CallType.FWD && full(lsn.callpair(cidx, 2)) == globalEidx
                        isForwardingTarget = true;
                        break;
                    end
                end

                if ~isCalledSync && ~isForwardingTarget
                    % Entry is only called asynchronously, skip reply-entry
                    continue;
                end

                if isempty(curEntry.replyActivity)
                    % the model was presumably loaded from a file without
                    % reply-activity specified by now lqns asks them so we
                    % need to calculate their location
                    racts = find(lsn.replygraph(:,eidx));
                    curEntry.replyActivity{1,end+1} = nodeHashMap{lsn.names{lsn.ashift+ racts}};
                end
                entryReplyElement = doc.createElement('reply-entry');
                taskActsElement.appendChild(entryReplyElement);
                entryReplyElement.setAttribute('name', nodeHashMap{curEntry.name});
                for r=1:length(curEntry.replyActivity)
                    entryReplyActElement = doc.createElement('reply-activity');
                    entryReplyElement.appendChild(entryReplyActElement);
                    entryReplyActElement.setAttribute('name', nodeHashMap{curEntry.replyActivity{r}});
                end
            end
        end
    end
end

%write the content into xml file
transformerFactory = TransformerFactory.newInstance();
transformer = transformerFactory.newTransformer();
transformer.setOutputProperty(OutputKeys.INDENT, 'yes');
transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
transformer.setOutputProperty(OutputKeys.VERSION, "1.0");
transformer.setOutputProperty(OutputKeys.ENCODING, "UTF-8");
transformer.setOutputProperty(OutputKeys.STANDALONE, "no");
source = DOMSource(doc);
if isempty(fileparts(filename))
    filename=[lineRootFolder,filesep,'workspace',filesep,filename];
end

%if self.options.verbose
%    line_printf('\nLQN model: %s\n', filename);
%end
result = StreamResult(File(filename));
transformer.transform(source, result);
end


function s = repl_to_lqnx(id)
% S = REPL_TO_LQNX(ID)
% ReplacementStrategy id -> the wire enum name, spelled exactly as the JSON
% interchange spells it (linemodel_save's repl_to_str). ONE spelling table for
% both transports: a second one would let a model round-trip through JSON and
% not through XML, or vice versa.
if id == ReplacementStrategy.LRU,       s = 'LRU';
elseif id == ReplacementStrategy.FIFO,  s = 'FIFO';
elseif id == ReplacementStrategy.RR,    s = 'RR';
elseif id == ReplacementStrategy.SFIFO, s = 'SFIFO';
elseif id == ReplacementStrategy.HLRU,  s = 'HLRU';
elseif id == ReplacementStrategy.CLIMB, s = 'CLIMB';
elseif id == ReplacementStrategy.QLRU,  s = 'QLRU';
else
    line_error(mfilename, sprintf('Unrecognized replacement strategy id %d.', id));
end
end


function s = lqnx_sched(schedText)
% S = LQNX_SCHED(SCHEDTEXT)
% The scheduling attribute lqns accepts for a discipline LINE holds by name.
%
% FCFSPRPRIO is LINE's reading of the LQN 'pri' discipline (LQIO SCHEDULE_PPR,
% preemptive priority resume), and 'fcfsprprio' is not a spelling lqns knows,
% so a round-tripped file would be rejected. Same mapping as the JAR
% lqnSchedText, the python sched_to_text and the C++ sched_to_lqnx.
if strcmpi(schedText, SchedStrategy.toText(SchedStrategy.FCFSPRPRIO))
    s = 'pri';
else
    s = schedText;
end
end


function s = lqnx_num(v)
% S = LQNX_NUM(V)
% A value attribute, at the shortest spelling that reads back as V.
%
% NOT num2str, which keeps five SIGNIFICANT digits: a demand or a call mean of
% 1/3 left as 0.33333, so three such calls summed to 0.99999 and a model
% round-tripped through .lqnx was a different model in the fifth digit. The
% file is a wire format for a solver, not a display, so it carries the value.
% Same rule as the JAR (Double.toString), python (repr) and the C++ writer's
% lqnx_num.
if isinf(v)
    if v > 0, s = 'Inf'; else, s = '-Inf'; end
    return
end
if isnan(v)
    s = 'NaN';
    return
end
for p = 15:17
    s = sprintf('%.*g', p, v);
    if sscanf(s, '%g') == v
        return
    end
end
end


function s = callGroupStrategyName(strategy)
% S = CALLGROUPSTRATEGYNAME(STRATEGY)
% RoutingStrategy id -> the wire enum name, spelled as the JSON interchange
% spells it. Only the two strategies a call group can be built with are
% named: WRROBIN would need per-target weights the group API does not take,
% and the remaining strategies are not dispatch policies at all, so an
% unnamed one is an error rather than a silent PROB.
if strategy == RoutingStrategy.RROBIN,  s = 'RROBIN';
elseif strategy == RoutingStrategy.JSQ, s = 'JSQ';
else
    line_error(mfilename, sprintf(['Call groups carry RROBIN or JSQ; routing ', ...
        'strategy %d cannot be written to .lqnx.'], strategy));
end
end


function popElement = popularity_element(doc, popularity, cardinality, entryName)
% POPELEMENT = POPULARITY_ELEMENT(DOC, POPULARITY, CARDINALITY, ENTRYNAME)
% <access-popularity name="CLASS"> with <parameter value="..."/> children in
% CONSTRUCTOR ORDER, the same ordering dist2json uses, so a popularity that
% round-trips through JSON round-trips here.
%
% A flat parameter list cannot by itself say where a DiscreteSampler's
% probability vector p ends and its support x begins, so the split is taken on
% the cardinality ALREADY DECLARED on the parent <item-entry>: n parameters
% means p over the default support 1..n, 2n means p followed by x. The default
% support is therefore OMITTED, being reconstructible.
%
% A class this encoding cannot carry is refused BY NAME. Writing a nameless
% element without its parameters is what this whole extension exists to stop.
popElement = doc.createElement('access-popularity');
cls = class(popularity);
popElement.setAttribute('name', cls);
switch cls
    case 'DiscreteSampler'
        p = popularity.getParam(1).paramValue(:)';
        x = popularity.getParam(2).paramValue(:)';
        vals = p;
        if numel(x) ~= numel(p) || any(abs(x - (1:numel(p))) > GlobalConstants.FineTol)
            vals = [p, x];
        end
    case 'Zipf'
        % Zipf(s, n): the shape and the item count, which is what the
        % constructor takes; p and x are derived from them.
        vals = [popularity.getParam(3).paramValue, popularity.getParam(4).paramValue];
    otherwise
        line_error(mfilename, sprintf(['the LQN XML interchange cannot encode a popularity ' ...
            'distribution of class ''%s'' on entry %s: <access-popularity> carries a flat ' ...
            'constructor-order parameter list, which this class does not have. Use ' ...
            'DiscreteSampler or Zipf, or extend popularity_element/popularity_from_element ' ...
            'in both writeXML.m and parseXML.m together.'], cls, entryName));
end
for v = 1:numel(vals)
    paramElement = doc.createElement('parameter');
    popElement.appendChild(paramElement);
    paramElement.setAttribute('value', lqnx_num(vals(v)));
end
end
