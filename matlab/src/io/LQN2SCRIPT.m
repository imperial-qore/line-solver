function LQN2SCRIPT(lqnmodel, modelName, fid)
% LQN2SCRIPT(MODEL, MODELNAME, FID)

% Copyright (c) 2012-2023, Imperial College London
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
fprintf(fid,'model = LayeredNetwork(''%s'');\n',modelName);
fprintf(fid,'\n');
%% host processors
for h=1:sn.nhosts
    fprintf(fid, 'P{%d} = Processor(model, ''%s'', %d, %s);\n', h, sn.names{h}, sn.mult(h), strrep(SchedStrategy.toFeature(sn.sched{h}),'_','.'));
    if sn.repl(h)~=1
        fprintf(fid, 'P{%d}.setReplication(%d);\n', h, sn.repl(h));
    end
end
fprintf(fid,'\n');
%% tasks
for t=1:sn.ntasks
    tidx = sn.tshift+t;
    fprintf(fid, 'T{%d} = Task(model, ''%s'', %d, %s).on(P{%d})', t, sn.names{tidx}, sn.mult(tidx), strrep(SchedStrategy.toFeature(sn.sched{tidx}),'_','.'),sn.parent(tidx));
    if sn.repl(tidx)~=1
        fprintf(fid, 'T{%d}.setReplication', t, sn.repl(tidx));
    end
    if ~isempty(sn.think{tidx})
        switch sn.think{tidx}.name
            case 'Immediate'
                fprintf(fid, '.setThinkTime(Immediate())');
            case 'Exp'
                fprintf(fid, '.setThinkTime(Exp.fitMean(%g))', sn.think{tidx}.getMean);
            case {'Erlang','HyperExp', 'Coxian', 'APH'}
                fprintf(fid, '.setThinkTime(%s.fitMeanAndSCV(%g,%g))', sn.think{tidx}.name, sn.think{tidx}.getMean, sn.think{tidx}.getSCV);
            otherwise
                line_error(mfilename,sprintf('LQN2SCRIPT does not support the %d distribution yet.',sn.think{tidx}.name));
        end
    end
    fprintf(fid,';\n');
end
fprintf(fid,'\n');
%% entries
for e=1:sn.nentries
    eidx = sn.eshift+e;
    fprintf(fid, 'E{%d} = Entry(model, ''%s'').on(T{%d});\n', e, sn.names{eidx},sn.parent(eidx)-sn.tshift);
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
        boundToStr = sprintf('.boundTo(E{%d})',boundTo);
    end
    if sn.sched{tidx} ~= SchedStrategy.ID_REF % ref tasks don't reply
        repliesToStr = '';
        repliesTo = find(sn.replygraph(aidx,:))-sn.eshift; % index of entry
        if ~isempty(repliesTo)
            if ~sn.isref(sn.parent(sn.eshift+repliesTo))
                repliesToStr = sprintf('.repliesTo(E{%d})',repliesTo);
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
                    callStr = sprintf('%s.synchCall(E{%d},%g)',callStr,calls(c)-sn.eshift,sn.callproc{c}.getMean);
                case CallType.ID_ASYNC
                    callStr = sprintf('%s.asynchCall(E{%d},%g)',callStr,calls(c)-sn.eshift,sn.callproc{c}.getMean);
            end
        end
    end
    switch sn.hostdem{aidx}.name
        case 'Immediate'
            fprintf(fid, 'A{%d} = Activity(model, ''%s'', Immediate()).on(T{%d})%s%s%s;\n', a, sn.names{aidx}, onTask, boundToStr, callStr, repliesToStr);
        case 'Exp'
            fprintf(fid, 'A{%d} = Activity(model, ''%s'', Exp.fitMean(%g)).on(T{%d})%s%s%s;\n', a, sn.names{aidx},sn.hostdem{aidx}.getMean, onTask, boundToStr, callStr, repliesToStr);
        case {'Erlang','HyperExp','Coxian','APH'}
            fprintf(fid, 'A{%d} = Activity(model, ''%s'', %s.fitMeanAndSCV(%g,%g)).on(T{%d})%s%s%s;\n', a, sn.names{aidx},sn.hostdem{aidx}.name,sn.hostdem{aidx}.getMean,sn.hostdem{aidx}.getSCV, onTask, boundToStr, callStr, repliesToStr);
        otherwise
            line_error(mfilename,sprintf('LQN2SCRIPT does not support the %d distribution yet.',sn.hostdem{aidx}.name));
    end
end
fprintf(fid,'\n');
%% think times
for h=1:sn.nhosts
    if ~isempty(sn.think{h})
        switch sn.think{h}.name
            case 'Immediate'
                fprintf(fid, 'P{%d}.setThinkTime(Immediate());\n', h);
            case 'Exp'
                fprintf(fid, 'P{%d}.setThinkTime(Exp.fitMean(%g));\n', h, sn.think{h}.getMean);
            case {'Erlang','HyperExp','Coxian','APH'}
                fprintf(fid, 'P{%d}.setThinkTime(%s.fitMeanAndSCV(%g,%g));\n', h, sn.think{h}.name, sn.think{h}.getMean, sn.think{h}.getSCV);
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
                fprintf(fid, 'T{%d}.addPrecedence(ActivityPrecedence.Serial(A{%d}, A{%d}));\n', tidx-sn.tshift, aidx-sn.ashift, bidx-sn.ashift);
            end
        end
    end
end

%% Loop precedences (POST_LOOP)
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
                        precActs = sprintf('A{%d}', bidx-sn.ashift);
                    else
                        precActs = sprintf('%s, A{%d}',precActs, bidx-sn.ashift);
                    end
                    if aidx ~= bidx % loop end reached
                        counts = 1/sn.graph(aidx,bidx);
                        fprintf(fid, 'T{%d}.addPrecedence(ActivityPrecedence.Loop(A{%d}, {%s}, %g));\n', tidx-sn.tshift, precMarker, precActs, counts);
                        precMarker = 0;
                    end
                end
            end
        end
    end
end

%% OrFork precedences (POST_OR)
precMarker = 0;
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    % for all successors
    for bidx=find(sn.graph(aidx,:))
        if bidx > sn.ashift % ignore precedence between entries and activities
            % Or pattern (POST_OR)
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_OR
                if precMarker == 0 % start a new orjoin
                    precMarker = aidx-sn.ashift;
                    precActs = sprintf('A{%d}', bidx-sn.ashift);
                    probs = sprintf('%g',full(sn.graph(aidx,bidx)));
                else
                    precActs = sprintf('%s, A{%d}', precActs, bidx-sn.ashift);
                    probs = sprintf('%s,%g',probs,full(sn.graph(aidx,bidx)));
                end
            end
        end
    end    
    if precMarker > 0
        fprintf(fid, 'T{%d}.addPrecedence(ActivityPrecedence.OrFork(A{%d},{%s},[%s]));\n', tidx-sn.tshift, precMarker, precActs, probs);
        precMarker = 0;
    end
end

%% AndFork precedences (POST_AND)
precMarker = 0;
for ai = 1:sn.nacts
    aidx = sn.ashift + ai;
    tidx = sn.parent(aidx);
    % for all successors
    for bidx=find(sn.graph(aidx,:))
        if bidx > sn.ashift % ignore precedence between entries and activities
            % Or pattern (POST_AND)
            if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_AND
                if precMarker == 0 % start a new orjoin
                    precMarker = aidx-sn.ashift;
                    precActs = sprintf('A{%d}', bidx-sn.ashift);
                else
                    precActs = sprintf('%s, A{%d}', precActs, bidx-sn.ashift);
                end
            end
        end
    end    
    if precMarker > 0
        fprintf(fid, 'T{%d}.addPrecedence(ActivityPrecedence.AndFork(A{%d},{%s}));\n', tidx-sn.tshift, precMarker, precActs);
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
                    precActs = sprintf('A{%d}', aidx-sn.ashift);
                else
                    precActs = sprintf('%s, A{%d}',precActs, aidx-sn.ashift);
                end
            end
        end
    end    
    if precMarker > 0
        fprintf(fid, 'T{%d}.addPrecedence(ActivityPrecedence.OrJoin({%s}, A{%d}));\n', tidx-sn.tshift, precActs, precMarker);
        precMarker = 0;
    end
end

%% OrJoin precedences (PRE_AND)
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
                    precActs = sprintf('A{%d}', aidx-sn.ashift);
                else
                    precActs = sprintf('%s, A{%d}',precActs, aidx-sn.ashift);
                end
            end
        end
    end    
    if precMarker > 0
        fprintf(fid, 'T{%d}.addPrecedence(ActivityPrecedence.AndJoin({%s}, A{%d}));\n', tidx-sn.tshift, precActs, precMarker);
        precMarker = 0;
    end
end

if ischar(fid)
    fclose(fid);
end
end
