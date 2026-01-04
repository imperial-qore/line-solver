function ind = getNodeIndex(self, name)
% NODEINDEX = GETNODEINDEX(NAME)

if isa(name,'Node')
    %node = name;
    %name = node.getName();
    ind = name.index;
    return
end
ind = find(cellfun(@(c) strcmp(c,name),self.getNodeNames));
end
