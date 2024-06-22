function LQN2JAVA(lqnmodel, modelName, fid)
% LQN2JAVA(MODEL, MODELNAME, FID)

% Copyright (c) 2012-2024, Imperial College London
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
fprintf(fid,'import jline.lang.distributions.*;\n');
fprintf(fid,'import jline.solvers.ln.SolverLN;\n\n');
fprintf(fid,'public class TestSolver%s {\n\n',upper(strrep(modelName,' ','')));
fprintf(fid,'\tpublic static void main(String[] args) throws Exception{\n\n');
fprintf(fid,'\tLayeredNetwork model = new LayeredNetwork("%s");\n',modelName);
fprintf(fid,'\n');
%% host processors
for h=1:sn.nhosts
    if isinf(sn.mult(h))
        fprintf(fid, '\tProcessor P%d = new Processor(model, "%s", Integer.MAX_VALUE, %s);\n', h, sn.names{h}, strrep(SchedStrategy.toFeature(sn.sched{h}),'_','.'));
    else
        fprintf(fid, '\tProcessor P%d = new Processor(model, "%s", %d, %s);\n', h, sn.names{h}, sn.mult(h), strrep(SchedStrategy.toFeature(sn.sched{h}),'_','.'));
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
        fprintf(fid, '\tTask T%d = new Task(model, "%s", Integer.MAX_VALUE, %s).on(P%d);\n', t, sn.names{tidx}, strrep(SchedStrategy.toFeature(sn.sched{tidx}),'_','.'), sn.parent(tidx));
    else
        fprintf(fid, '\tTask T%d = new Task(model, "%s", %d, %s).on(P%d);\n', t, sn.names{tidx}, sn.mult(tidx), strrep(SchedStrategy.toFeature(sn.sched{tidx}),'_','.'), sn.parent(tidx));
    end
    if sn.repl(tidx)~=1
        fprintf(fid, '\tT%d.setReplication(%d);\n', t, sn.repl(tidx));
    end
    if ~isempty(sn.think{tidx})
        switch sn.think{tidx}.name
            case 'Immediate'
                fprintf(fid, '\tT%d.setThinkTime(new Immediate());\n',t);
            case 'Exp'
                fprintf(fid, '\tT%d.setThinkTime(new Exp(%g));\n',t, 1/sn.think{tidx}.getMean);
            case {'Erlang','HyperExp', 'Coxian', 'APH'}
                fprintf(fid, '\tT%d.setThinkTime(%s.fitMeanAndSCV(%g,%g));\n', t, sn.think{tidx}.name, sn.think{tidx}.getMean, sn.think{tidx}.getSCV);
            otherwise
                line_error(mfilename,sprintf('LQN2SCRIPT does not support the %d distribution yet.',sn.think{tidx}.name));
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
    if sn.sched{tidx} ~= SchedStrategy.ID_REF % ref tasks don't reply
        repliesToStr = '';
        repliesTo = find(sn.replygraph(aidx,:))-sn.eshift; % index of entry
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
                case CallType.ID_SYNC
                    callStr = sprintf('%s.synchCall(E%d,%g)',callStr,calls(c)-sn.eshift,sn.callproc{c}.getMean);
                case CallType.ID_ASYNC
                    callStr = sprintf('%s.asynchCall(E%d,%g)',callStr,calls(c)-sn.eshift,sn.callproc{c}.getMean);
            end
        end
    end
    switch sn.hostdem{aidx}.name
        case 'Immediate'
            fprintf(fid, '\tActivity A%d = new Activity(model, "%s", new Immediate()).on(T%d);', a, sn.names{aidx}, onTask);
        case 'Exp'
            fprintf(fid, '\tActivity A%d = new Activity(model, "%s", new Exp(%g)).on(T%d);', a, sn.names{aidx},1/sn.hostdem{aidx}.getMean, onTask);
        case {'Erlang','HyperExp','Coxian','APH'}
            fprintf(fid, '\tActivity A%d = new Activity(model, "%s", %s.fitMeanAndSCV(%g,%g)).on(T%d);', a, sn.names{aidx},sn.hostdem{aidx}.name,sn.hostdem{aidx}.getMean,sn.hostdem{aidx}.getSCV, onTask);
        otherwise
            line_error(mfilename,sprintf('LQN2SCRIPT does not support the %d distribution yet.',sn.hostdem{aidx}.name));
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
    if ~isempty(sn.think{h})
        switch sn.think{h}.name
            case 'Immediate'
                fprintf(fid, '\tP%d.setThinkTime(Immediate());\n', h);
            case 'Exp'
                fprintf(fid, '\tP%d.setThinkTime(Exp(%g));\n', h, sn.think{h}.getMean);
            case {'Erlang','HyperExp','Coxian','APH'}
                fprintf(fid, '\tP%d.setThinkTime(%s.fitMeanAndSCV(%g,%g));\n', h, sn.think{h}.name, sn.think{h}.getMean, sn.think{h}.getSCV);
            otherwise
                line_error(mfilename,sprintf('LQN2SCRIPT does not support the %d distribution yet.',sn.think{h}.name));
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
            if full(sn.actpretype(aidx)) == ActivityPrecedenceType.ID_PRE_SEQ && full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_SEQ
                fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.Sequence("%s", "%s"));\n', tidx-sn.tshift, sn.names{aidx}, sn.names{bidx});
            end
        end
    end
end

%% Loop precedences (POST_LOOP)
hasPreActs = false;
hasPostActs = false;
precMarker = 0;
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    % for all successors
    for bidx=find(sn.graph(aidx,:))
        if bidx > sn.ashift % ignore precedence between entries and activities
            % Loop pattern (POST_LOOP)
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_LOOP
                if precMarker == 0 % start a new loop
                    precMarker = aidx-sn.ashift;
                    precActs = '';
                else
                    if isempty(precActs)
                        precActs = sprintf('\tprecActs.add("%s");\n', sn.names{bidx});
                    else
                        precActs = sprintf('%s\tprecActs.add("%s");\n', precActs, sn.names{bidx});
                    end
                    if aidx ~= bidx % loop end reached
                        counts = 1/sn.graph(aidx,bidx);
                        fprintf(fid, '\n\t// Loop Activity Precedence \n');
                        if ~hasPreActs
                            fprintf(fid, '\tArrayList<String> precActs = new ArrayList<String>();\n');
                            hasPreActs = true;
                        else
                            fprintf(fid, '\tprecActs.clear();\n');
                        end
                        fprintf(fid, precActs);
                        fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.Loop("%s", precActs, new Matrix(%g)));\n', tidx-sn.tshift, sn.names{precMarker+sn.ashift}, counts);
                        precMarker = 0;
                    end
                end
            end
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
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_OR
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
            fprintf(fid, '\tprecActs.clear();\n');
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
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_AND
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
            fprintf(fid, '\t postActs.clear();\n');
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
            if full(sn.actpretype(aidx)) == ActivityPrecedenceType.ID_PRE_OR
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
            fprintf(fid, '\tprecActs.clear();\n');
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
            if full(sn.actpretype(aidx)) == ActivityPrecedenceType.ID_PRE_AND
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
            fprintf(fid, '\tprecActs.clear();\n');
        end
        fprintf(fid, precActs);
        fprintf(fid, '\tT%d.addPrecedence(ActivityPrecedence.AndJoin(precActs, "%s"));\n', tidx-sn.tshift, sn.names{precMarker+sn.ashift});
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
