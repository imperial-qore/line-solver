function clone = copyElement(self)
% CLONE = COPYELEMENT()

% Make a shallow copy of all properties
clone = copyElement@Copyable(self);
% Make a deep copy of each handle
for i=1:length(self.classes)
    clone.classes{i} = self.classes{i}.copy;
end
% Make a deep copy of each handle
for i=1:length(self.nodes)
    clone.nodes{i} = self.nodes{i}.copy;
    if isa(clone.nodes{i},'Station')
        clone.stations{i} = clone.nodes{i};
    end
    clone.connections = self.connections;
end
end
