function [taggedModel, taggedJob] = tagChain(model, chain, jobclass, suffix)
% the tagged job will be removed from the initial
% population of JOBCLASS
if nargin<4 || isempty(suffix)
    suffix = '.tagged';
end
if nargin<3
    jobclass = chain.classes{1};
end
I = model.getNumberOfNodes;
R = model.getNumberOfClasses;
taggedModel = model.copy;

% we don't use rtNodesByClass because it contains the
% fictitious class switching nodes
Plinked = taggedModel.getLinkedRoutingMatrix;
if ~iscell(Plinked) || isempty(Plinked)
    line_error(mfilename, 'getCdfRespT requires the original model to be linked with a routing matrix defined as a cell array P{r,s} for every class pair (r,s).');
end

taggedModel.resetNetwork; % resets cs Nodes as well
taggedModel.reset(true);

chainIndexes = cell2mat(chain.index);
for r=chainIndexes
    % create a tagged class
    taggedModel.classes{end+1,1} = taggedModel.classes{r}.copy;
    taggedModel.classes{end,1}.index=length(taggedModel.classes);
    taggedModel.classes{end,1}.name=[taggedModel.classes{r,1}.name,suffix];
    if r==jobclass.index
        taggedModel.classes{end}.population = 1;
    else
        taggedModel.classes{end}.population = 0;
    end

    % clone station sections for tagged class
    for m=1:length(taggedModel.nodes)
        taggedModel.stations{m}.output.outputStrategy{end+1} = taggedModel.stations{m}.output.outputStrategy{r};
    end

    for m=1:length(taggedModel.stations)
        if model.stations{m}.server.serviceProcess{r}{end}.isDisabled
            taggedModel.stations{m}.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{end}{end}.copy;
            taggedModel.stations{m}.server.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r};
            taggedModel.stations{m}.server.serviceProcess{end}{end}=taggedModel.stations{m}.server.serviceProcess{r}{end}.copy;
            taggedModel.stations{m}.schedStrategyPar(end+1) = 0;
            taggedModel.stations{m}.dropRule(1,end+1) = -1;
            taggedModel.stations{m}.classCap(1,r) = 0;
            taggedModel.stations{m}.classCap(1,end+1) = 0;
        else
            taggedModel.stations{m}.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r}{end}.copy;
            taggedModel.stations{m}.server.serviceProcess{end+1} = taggedModel.stations{m}.server.serviceProcess{r};
            taggedModel.stations{m}.server.serviceProcess{end}{end}=taggedModel.stations{m}.server.serviceProcess{r}{end}.copy;
            taggedModel.stations{m}.schedStrategyPar(end+1) = taggedModel.stations{m}.schedStrategyPar(r);
            taggedModel.stations{m}.classCap(1,end+1) = 1;
            taggedModel.stations{m}.dropRule(1,end+1) = -1;
            taggedModel.stations{m}.classCap(1,r) = taggedModel.stations{m}.classCap(r) - 1;
        end
    end
end

taggedModel.classes{jobclass.index,1}.population = taggedModel.classes{jobclass.index}.population - 1;

for ir=1:length(chainIndexes)
    r = chainIndexes(ir);
    for is=1:length(chainIndexes)
        s = chainIndexes(is);
        Plinked{R+ir,R+is} = Plinked{r,s};
    end
end
Rp = taggedModel.getNumberOfClasses;
for r=1:Rp
    for s=1:Rp
        if isempty(Plinked{r,s})
            Plinked{r,s} = zeros(I);
        end
    end
end
taggedModel.sn = [];
taggedModel.link(Plinked);
taggedModel.reset(true);
taggedModel.refreshStruct(true);
taggedModel.initDefault;
tchains = taggedModel.getChains;
taggedJob = tchains{end};
end