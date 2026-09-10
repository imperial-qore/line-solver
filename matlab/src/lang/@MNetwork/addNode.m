function isReplacement = addNode(self, node)
% ADDNODE(NODE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%
hasSameName = cellfun(@(x) strcmp(x.name,node.name), {self.nodes{1:end}});
if any(hasSameName)
    if self.allowReplace
        oldNode = self.nodes{find(hasSameName)}; %#ok<FNDSB>
        node.index = oldNode.index;
        self.nodes{oldNode.index,1} = node;
        if isa(node,'Station')
            self.stations{oldNode.stationIndex,1} = node;
            node.stationIndex = oldNode.stationIndex;
        end
        self.setUsedLangFeature(class(node)); % station type
        self.setUsedLangFeature(class(node.input)); % station type
        self.setUsedLangFeature(class(node.server)); % station type
        self.setUsedLangFeature(class(node.output)); % station type
        isReplacement = true;
    else
        line_error(mfilename,'A node with an identical name already exists.');
    end
    return;
end
self.nodes{end+1,1} = node;
node.index = length(self.nodes);
if isa(node,'Station')
    self.stations{end+1,1} = node;
    node.stationIndex = length(self.stations);
end
self.setUsedLangFeature(class(node)); % station type
self.setUsedLangFeature(class(node.input)); % station type
self.setUsedLangFeature(class(node.server)); % station type
self.setUsedLangFeature(class(node.output)); % station type
isReplacement = false;
end
