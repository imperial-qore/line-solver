function [wf, actIdxOf, callIdxOf, execs, callexecs] = lqn_entry_workflow(model, lqn, eidx, withCalls)
% [WF, ACTIDXOF, CALLIDXOF, EXECS, CALLEXECS] = LQN_ENTRY_WORKFLOW(MODEL, LQN, EIDX, WITHCALLS)
%
% Activity graph of LQN entry EIDX as a Workflow, so that the entry service
% time can be composed exactly into a phase-type law by the series-parallel
% reduction of Workflow.toPH.
%
% MODEL is the LayeredNetwork the struct LQN was obtained from; the activity
% precedences are read from the task object rather than reconstructed from
% LQN.GRAPH, whose loop back-edges carry probabilities and not counts.
%
% WITHCALLS true expands every synchronous call of an activity into a leaf of
% its own, placed in series after the activity, so that the call response law
% and the host demand law stay separable across iterations. WITHCALLS false
% keeps only the host demands, which is the processor-demand law of the entry:
% the host is released while a call is outstanding.
%
% Returns:
%   wf        - Workflow whose leaves are the activities (and calls) of EIDX
%   actIdxOf  - (nidx,1) workflow activity index of each LQN activity, 0 if absent
%   callIdxOf - (ncalls,1) workflow activity index of each call, 0 if absent
%   execs     - (nidx,1) expected executions of each LQN activity per entry
%               invocation
%   callexecs - (ncalls,1) expected executions of each call per entry invocation,
%               not counting the mean number of calls per execution
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 4
    withCalls = true;
end

tidx = lqn.parent(eidx);
t = tidx - lqn.tshift;
if t < 1 || t > numel(model.tasks)
    line_error(mfilename, sprintf('Entry %s has no task in the layered model.', lqn.hashnames{eidx}));
end
task = model.tasks{t};

acts = lqn.actsof{eidx};
acts = acts(:)';
if isempty(acts)
    line_error(mfilename, sprintf('Entry %s binds no activity.', lqn.hashnames{eidx}));
end

wf = Workflow([lqn.names{eidx}, '.Workflow']);

actIdxOf = zeros(lqn.nidx, 1);
callIdxOf = zeros(lqn.ncalls, 1);
headName = cell(lqn.nidx, 1);
tailName = cell(lqn.nidx, 1);
% name of an activity of this entry -> its absolute index, used to rewrite the
% precedences of the task onto the expanded chains
actOfName = cell(0, 2);

for aidx = acts
    nm = lqn.names{aidx};
    a = wf.addActivity(nm, lqn.hostdem{aidx});
    actIdxOf(aidx) = a.index;
    headName{aidx} = nm;
    tailName{aidx} = nm;
    actOfName(end+1, 1:2) = {nm, aidx}; %#ok<AGROW>
end

if withCalls
    for aidx = acts
        chain = {headName{aidx}};
        for cidx = lqn.callsof{aidx}
            if lqn.calltype(cidx) ~= CallType.SYNC
                continue % an asynchronous call blocks the caller for no time
            end
            cnm = lqn.callhashnames{cidx};
            c = wf.addActivity(cnm, Immediate.getInstance());
            callIdxOf(cidx) = c.index;
            chain{end+1} = cnm; %#ok<AGROW>
        end
        if numel(chain) > 1
            wf.addPrecedence(Workflow.Serial(chain{:}));
            tailName{aidx} = chain{end};
        end
    end
end

% Precedences of the task, restricted to the activities of this entry and
% rewritten so that a predecessor is entered at its head and left at its tail
for p = 1:length(task.precedences)
    prec = task.precedences(p);
    preIdx = resolveNames(prec.preActs, actOfName);
    postIdx = resolveNames(prec.postActs, actOfName);
    if any(preIdx == 0) || any(postIdx == 0)
        continue % the precedence belongs to another entry of the same task
    end
    preNames = cell(1, numel(preIdx));
    for k = 1:numel(preIdx)
        preNames{k} = tailName{preIdx(k)};
    end
    postNames = cell(1, numel(postIdx));
    for k = 1:numel(postIdx)
        postNames{k} = headName{postIdx(k)};
    end
    wf.addPrecedence(ActivityPrecedence(preNames, postNames, ...
        prec.preType, prec.postType, prec.preParams, prec.postParams));
end

[isValid, msg] = wf.validate();
if ~isValid
    line_error(mfilename, sprintf('Entry %s cannot be composed into a phase-type law: %s', ...
        lqn.hashnames{eidx}, msg));
end

tree = wf.getSPTree();
if isempty(tree)
    line_error(mfilename, sprintf(['Entry %s has a precedence graph that is not series-parallel, ' ...
        'so its activity graph has no exact phase-type reduction. Use method=''default''.'], ...
        lqn.hashnames{eidx}));
end

execs = zeros(lqn.nidx, 1);
for aidx = acts
    execs(aidx) = tree.execs(tree.leafOf(actIdxOf(aidx)));
end
callexecs = zeros(lqn.ncalls, 1);
for cidx = find(callIdxOf > 0)'
    callexecs(cidx) = tree.execs(tree.leafOf(callIdxOf(cidx)));
end
end

function idx = resolveNames(names, actOfName)
% IDX = RESOLVENAMES(NAMES, ACTOFNAME) maps activity names to absolute indices,
% returning 0 for a name outside this entry
idx = zeros(1, numel(names));
for k = 1:numel(names)
    nm = names{k};
    if isa(nm, 'Activity')
        nm = nm.getName();
    end
    for r = 1:size(actOfName, 1)
        if strcmp(actOfName{r,1}, nm)
            idx(k) = actOfName{r,2};
            break
        end
    end
end
end
