function LQN2JAVA(lqnmodel, modelName, fid)
% LQN2JAVA(MODEL, MODELNAME, FID)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2%~exist('modelName','var')
    modelName='myLayeredModel';
end
if nargin<3%~exist('fid','var')
    fid=1;
end
if ischar(fid)
    fid = fopen(fid,'w');
end
sn = lqnmodel.getStruct;

%% initialization
fprintf(fid,'package jline.examples;\n\n');
fprintf(fid,'import jline.lang.*;\n');
fprintf(fid,'import jline.lang.constant.*;\n');
fprintf(fid,'import jline.lang.processes.*;\n');
fprintf(fid,'import jline.solvers.ln.SolverLN;\n\n');
fprintf(fid,'public class TestSolver%s {\n\n',upper(strrep(modelName,' ','')));
fprintf(fid,'\tpublic static void main(String[] args) throws Exception{\n\n');
fprintf(fid,'\tLayeredNetwork model = new LayeredNetwork("%s");\n',modelName);
fprintf(fid,'\n');
%% host processors
for h=1:sn.nhosts
    if isinf(sn.mult(h))
        fprintf(fid, '\tProcessor P%d = new Processor(model, "%s", Integer.MAX_VALUE, %s);\n', h, sn.names{h}, strrep(SchedStrategy.toFeature(sn.sched(h)),'_','.'));
    else
        fprintf(fid, '\tProcessor P%d = new Processor(model, "%s", %d, %s);\n', h, sn.names{h}, sn.mult(h), strrep(SchedStrategy.toFeature(sn.sched(h)),'_','.'));
    end
    if sn.repl(h)~=1
        fprintf(fid, 'P%d.setReplication(%d);\n', h, sn.repl(h));
    end
end
fprintf(fid,'\n');
%% tasks
for t=1:sn.ntasks
    tidx = sn.tshift+t;
    if isinf(sn.mult(tidx))
        fprintf(fid, '\tTask T%d = new Task(model, "%s", Integer.MAX_VALUE, %s).on(P%d);\n', t, sn.names{tidx}, strrep(SchedStrategy.toFeature(sn.sched(tidx)),'_','.'), sn.parent(tidx));
    else
        fprintf(fid, '\tTask T%d = new Task(model, "%s", %d, %s).on(P%d);\n', t, sn.names{tidx}, sn.mult(tidx), strrep(SchedStrategy.toFeature(sn.sched(tidx)),'_','.'), sn.parent(tidx));
    end
    if sn.repl(tidx)~=1
        fprintf(fid, '\tT%d.setReplication(%d);\n', t, sn.repl(tidx));
    end
    if ~isempty(sn.think{tidx}) && sn.think_type(tidx) ~= ProcessType.DISABLED
        switch sn.think_type(tidx)
            case ProcessType.IMMEDIATE
                fprintf(fid, '\tT%d.setThinkTime(new Immediate());\n',t);
            case ProcessType.EXP
                fprintf(fid, '\tT%d.setThinkTime(new Exp(%g));\n',t, 1/sn.think_mean(tidx));
            case {ProcessType.ERLANG, ProcessType.HYPEREXP, ProcessType.COXIAN, ProcessType.APH}
                fprintf(fid, '\tT%d.setThinkTime(%s.fitMeanAndSCV(%g,%g));\n', t, char(sn.think_type(tidx)), sn.think_mean(tidx), sn.think_scv(tidx));
            otherwise
                line_error(mfilename,sprintf('LQN2JAVA does not support the %s distribution yet.',char(sn.think_type(tidx))));
        end
    end
end
fprintf(fid,'\n');
%% entries
for e=1:sn.nentries
    eidx = sn.eshift+e;
    fprintf(fid, '\tEntry E%d = new Entry(model, "%s").on(T%d);\n', e, sn.names{eidx},sn.parent(eidx)-sn.tshift);
end
fprintf(fid,'\n');
%% activities
for a=1:sn.nacts
    aidx = sn.ashift+a;
    tidx = sn.parent(aidx);
    onTask = tidx-sn.tshift;
    boundTo = find(sn.graph((sn.eshift+1):(sn.eshift+sn.nentries),aidx));
    boundToStr = '';
    if ~isempty(boundTo)
        boundToStr = sprintf('.boundTo(E%d)',boundTo);
    end
    repliesToStr = '';    
    if sn.sched(tidx) ~= SchedStrategy.REF % ref tasks don't reply        
        repliesTo = find(sn.replygraph(a,:)); % index of entry
        if ~isempty(repliesTo)
            if ~sn.isref(sn.parent(sn.eshift+repliesTo))
                repliesToStr = sprintf('.repliesTo(E%d)',repliesTo);
            end
        end
    end
    callStr = '';
    if ~isempty(sn.callpair)
        cidxs = find(sn.callpair(:,1)==aidx);
        calls = sn.callpair(:,2);
        for c=cidxs(:)'
            switch sn.calltype(c)
                case CallType.SYNC
                    callStr = sprintf('%s.synchCall(E%d,%g)',callStr,calls(c)-sn.eshift,sn.callproc_mean(c));
                case CallType.ASYNC
                    callStr = sprintf('%s.asynchCall(E%d,%g)',callStr,calls(c)-sn.eshift,sn.callproc_mean(c));
            end
        end
    end
    switch sn.hostdem_type(aidx)
        case ProcessType.IMMEDIATE
            fprintf(fid, '\tActivity A%d = new Activity(model, "%s", new Immediate()).on(T%d);', a, sn.names{aidx}, onTask);
        case ProcessType.EXP
            fprintf(fid, '\tActivity A%d = new Activity(model, "%s", new Exp(%g)).on(T%d);', a, sn.names{aidx},1/sn.hostdem_mean(aidx), onTask);
        case {ProcessType.ERLANG, ProcessType.HYPEREXP, ProcessType.COXIAN, ProcessType.APH}
            fprintf(fid, '\tActivity A%d = new Activity(model, "%s", %s.fitMeanAndSCV(%g,%g)).on(T%d);', a, sn.names{aidx},char(sn.hostdem_type(aidx)),sn.hostdem_mean(aidx),sn.hostdem_scv(aidx), onTask);
        otherwise
            line_error(mfilename,sprintf('LQN2JAVA does not support the %s distribution yet.',char(sn.hostdem_type(aidx))));
    end
    if ~isempty(boundToStr)
        fprintf(fid, ' A%d%s;', a, boundToStr);
    end
    if ~isempty(callStr)
        fprintf(fid, ' A%d%s;', a, callStr);
    end
    if ~isempty(repliesToStr)
        fprintf(fid, ' A%d%s;', a, repliesToStr);
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
%% think times
for h=1:sn.nhosts
    if ~isempty(sn.think{h}) && sn.think_type(h) ~= ProcessType.DISABLED
        switch sn.think_type(h)
            case ProcessType.IMMEDIATE
                fprintf(fid, '\tP%d.setThinkTime(Immediate());\n', h);
            case ProcessType.EXP
                fprintf(fid, '\tP%d.setThinkTime(Exp(%g));\n', h, sn.think_mean(h));
            case {ProcessType.ERLANG, ProcessType.HYPEREXP, ProcessType.COXIAN, ProcessType.APH}
                fprintf(fid, '\tP%d.setThinkTime(%s.fitMeanAndSCV(%g,%g));\n', h, char(sn.think_type(h)), sn.think_mean(h), sn.think_scv(h));
            otherwise
                line_error(mfilename,sprintf('LQN2JAVA does not support the %s distribution yet.',char(sn.think_type(h))));
        end
    end
end

%% Sequential precedences
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    % for all successors
    for bidx=find(sn.graph(aidx,:))
        if bidx > sn.ashift % ignore precedence between entries and activities
            % Serial pattern (SEQ)
            if full(sn.actpretype(aidx)) == ActivityPrecedenceType.PRE_SEQ && full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_SEQ
                fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.Serial("%s", "%s"));\n', tidx-sn.tshift, sn.names{aidx}, sn.names{bidx});
            end
        end
    end
end

%% Loop precedences (POST_LOOP)
% Loop structure in sn.graph:
%   - Entry activity (preAct) has edge to first loop body activity
%   - Last loop body activity has back-edge to first loop body (weight = 1-1/counts)
%   - Last loop body activity has edge to end activity (weight = 1/counts)
%   - All loop body and end activities have POST_LOOP type
hasPreActs = false;
hasPostActs = false;
processedLoops = false(1, sn.nacts);
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    % Check if this activity starts a loop (has a successor with POST_LOOP type)
    % and hasn't been processed as part of another loop
    if processedLoops(ai)
        continue;
    end

    successors = find(sn.graph(aidx,:));
    for bidx = successors
        if bidx > sn.ashift && full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_LOOP
            % Skip if this loop body activity was already processed
            if processedLoops(bidx - sn.ashift)
                continue;
            end
            % Found start of a loop: aidx is the entry, bidx is first loop body activity
            loopStart = bidx;
            precMarker = ai;
            loopActNames = {};

            % Follow the chain of POST_LOOP activities
            curIdx = loopStart;
            while true
                loopActNames{end+1} = sn.names{curIdx}; %#ok<AGROW>
                processedLoops(curIdx - sn.ashift) = true;

                % Find successors of current activity
                curSuccessors = find(sn.graph(curIdx,:));
                curSuccessors = curSuccessors(curSuccessors > sn.ashift);

                % Check for loop termination: find the end activity
                % End activity has weight = 1/counts (not the back-edge weight)
                endIdx = 0;
                nextIdx = 0;
                for succIdx = curSuccessors
                    if full(sn.actposttype(succIdx)) == ActivityPrecedenceType.POST_LOOP
                        if succIdx == loopStart
                            % This is the back-edge, skip it
                            continue;
                        end
                        weight = full(sn.graph(curIdx, succIdx));
                        if weight > 0 && weight < 1
                            % This is the end activity (weight = 1/counts)
                            endIdx = succIdx;
                        else
                            % This is the next activity in the loop body (weight = 1.0)
                            nextIdx = succIdx;
                        end
                    end
                end

                if endIdx > 0
                    % Found end activity - calculate counts and output
                    weight = full(sn.graph(curIdx, endIdx));
                    if weight > 0
                        counts = 1/weight;
                    else
                        counts = 1; % Fallback to prevent division by zero
                    end
                    loopActNames{end+1} = sn.names{endIdx}; %#ok<AGROW>
                    processedLoops(endIdx - sn.ashift) = true;

                    % Build the precActs string
                    fprintf(fid, '\n\t// Loop Activity Precedence \n');
                    if ~hasPreActs
                        fprintf(fid, '\tArrayList<String> precActs = new ArrayList<String>();\n');
                        hasPreActs = true;
                    else
                        fprintf(fid, '\tprecActs = new ArrayList<String>();\n');
                    end
                    for k = 1:length(loopActNames)
                        fprintf(fid, '\tprecActs.add("%s");\n', loopActNames{k});
                    end
                    fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.Loop("%s", precActs, Matrix.singleton(%g)));\n', tidx-sn.tshift, sn.names{aidx}, counts);
                    break;
                elseif nextIdx > 0
                    % Continue to next activity in loop body
                    curIdx = nextIdx;
                else
                    % No more successors - shouldn't happen in valid loop
                    break;
                end
            end
            break; % Only process one loop starting from this activity
        end
    end
end

%% OrFork precedences (POST_OR)
hasProbs = false;
precMarker = 0;
precActs = '';
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    prob_ctr = 0;
    % for all successors
    for bidx=find(sn.graph(aidx,:))
        if bidx > sn.ashift % ignore precedence between entries and activities
            % Or pattern (POST_OR)
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_OR
                if precMarker == 0 % start a new orjoin
                    precActs = '';
                    precMarker = aidx-sn.ashift;
                    prob_ctr = prob_ctr + 1;
                    precActs = sprintf('\tprecActs.add("%s");\n', sn.names{bidx});
                    probString = sprintf('\tprobs.set(0,%d,%g);\n', prob_ctr-1, full(sn.graph(aidx,bidx)));
                else
                    prob_ctr = prob_ctr + 1;
                    precActs = sprintf('%s\tprecActs.add("%s");\n', precActs, sn.names{bidx});
                    probString = sprintf('%s\tprobs.set(0,%d,%g);\n', probString, prob_ctr-1, full(sn.graph(aidx,bidx)));
                end
            end
        end
    end


    if precMarker > 0
        fprintf(fid, '\n\t// OrFork Activity Precedence \n');
        if ~hasPreActs
            fprintf(fid, '\tArrayList<String> precActs = new ArrayList<String>();\n');
            hasPreActs = true;
        else
            fprintf(fid, '\tprecActs = new ArrayList<String>();\n');
        end
        if ~hasProbs
            fprintf(fid, '\tMatrix probs = new Matrix(1,%d);\n',prob_ctr);
            hasProbs = true;
        else
            fprintf(fid, '\tprobs = new Matrix(1,%d);\n',prob_ctr);
        end
        fprintf(fid, precActs);
        fprintf(fid, probString);
        fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.OrFork("%s", precActs, probs));\n', tidx-sn.tshift, sn.names{precMarker+sn.ashift});
        precMarker = 0;
    end
end

%% AndFork precedences (POST_AND)
precMarker = 0;
postActs = '';
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    % for all successors
    for bidx=find(sn.graph(aidx,:))
        if bidx > sn.ashift % ignore precedence between entries and activities
            % Or pattern (POST_AND)
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_AND
                if isempty(postActs)
                    postActs = sprintf('\tpostActs.add("%s");\n', sn.names{bidx});
                else
                    postActs = sprintf('%s\tpostActs.add("%s");\n', postActs, sn.names{bidx});
                end

                if precMarker == 0 % start a new orjoin
                    precMarker = aidx-sn.ashift;
                end
            end
        end
    end
    if precMarker > 0
        fprintf(fid, '\n\t// AndFork Activity Precedence \n');
        if ~hasPostActs
            fprintf(fid, '\tArrayList<String> postActs = new ArrayList<String>();\n');
            hasPostActs = true;
        else
            fprintf(fid, '\t postActs = new ArrayList<String>();\n');
        end

        fprintf(fid, postActs);
        fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.AndFork("%s", postActs));\n', tidx-sn.tshift, sn.names{precMarker+sn.ashift});
        precMarker = 0;
    end
end


%% OrJoin precedences (PRE_OR)
precMarker = 0;
for bi = sn.nacts:-1:1
    bidx = sn.ashift + bi;
    tidx = sn.parent(bidx);
    % for all predecessors
    for aidx=find(sn.graph(:,bidx))'
        if aidx > sn.ashift % ignore precedence between entries and activities
            % OrJoin pattern (PRE_OR)
            if full(sn.actpretype(aidx)) == ActivityPrecedenceType.PRE_OR
                if precMarker == 0 % start a new orjoin
                    precMarker = bidx-sn.ashift;
                    precActs = sprintf('\tprecActs.add("%s");\n', sn.names{aidx});
                else
                    precActs = sprintf('%s\tprecActs.add("%s");\n', precActs, sn.names{aidx});
                end
            end
        end
    end
    if precMarker > 0
        fprintf(fid, '\n\t// OrJoin Activity Precedence \n');
        if ~hasPreActs
            fprintf(fid, '\tArrayList<String> precActs = new ArrayList<String>();\n');
            hasPreActs = true;
        else
            fprintf(fid, '\tprecActs = new ArrayList<String>();\n');
        end

        fprintf(fid, precActs);
        fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.OrJoin(precActs, "%s"));\n', tidx-sn.tshift, sn.names{precMarker+sn.ashift});
        precMarker = 0;
    end
end

%% AndJoin precedences (PRE_AND)
precMarker = 0;
for bi = sn.nacts:-1:1
    bidx = sn.ashift + bi;
    tidx = sn.parent(bidx);
    % for all predecessors
    for aidx=find(sn.graph(:,bidx))'
        if aidx > sn.ashift % ignore precedence between entries and activities
            % OrJoin pattern (PRE_AND)
            if full(sn.actpretype(aidx)) == ActivityPrecedenceType.PRE_AND
                if precMarker == 0 % start a new orjoin
                    precMarker = bidx-sn.ashift;
                    precActs = sprintf('\tprecActs.add("%s");\n', sn.names{aidx});
                else
                    precActs = sprintf('%s\tprecActs.add("%s");\n', precActs, sn.names{aidx});
                end
            end
        end
    end
    if precMarker > 0
        fprintf(fid, '\n\t// AndJoin Activity Precedence \n');
        if ~hasPreActs
            fprintf(fid, '\tArrayList<String> precActs = new ArrayList<String>();\n');
            hasPreActs = true;
        else
            fprintf(fid, '\tprecActs = new ArrayList<String>();\n');
        end
        fprintf(fid, precActs);
        
        % Find quorum parameter from original precedence structure
        postActName = sn.names{precMarker+sn.ashift};
        localTaskIdx = tidx - sn.tshift;
        quorum = [];
        for ap = 1:length(lqnmodel.tasks{localTaskIdx}.precedences)
            precedence = lqnmodel.tasks{localTaskIdx}.precedences(ap);
            if precedence.preType == ActivityPrecedenceType.PRE_AND
                % Check if this precedence contains our post activity
                if any(strcmp(precedence.postActs, postActName))
                    quorum = precedence.preParams;
                    break;
                end
            end
        end
        
        if isempty(quorum)
            fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, "%s"));\n', localTaskIdx, postActName);
        else
            fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, "%s", %g));\n', localTaskIdx, postActName, quorum);
        end
        precMarker = 0;
    end
end

fprintf(fid, '\n\t// Model solution \n');
fprintf(fid,'\tSolverLN solver = new SolverLN(model);\n');
fprintf(fid,'\tsolver.getEnsembleAvg();\n');
fprintf(fid,'\t}\n}\n');

if ischar(fid)
    fclose(fid);
end
end
