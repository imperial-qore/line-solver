function [names,hostnames,tasknames,entrynames,actnames] = getNodeNames(self)
% NAME = GETNODENAMES(SELF)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
hostnames = cellfun(@(x)x.name,self.hosts,'UniformOutput',false);
tasknames = cellfun(@(x)x.name,self.tasks,'UniformOutput',false);
entrynames = cellfun(@(x)x.name,self.entries,'UniformOutput',false);
actnames = cellfun(@(x)x.name,self.activities,'UniformOutput',false);
names = [hostnames(:)',tasknames(:)',entrynames(:)',actnames(:)'];
end
