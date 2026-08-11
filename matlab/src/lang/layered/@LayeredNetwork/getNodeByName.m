function node = getNodeByName(self,name)
% NODE = GETNODEBYNAME(SELF,NAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
[~,hostnames,tasknames,entrynames,actnames] = getNodeNames(self);

% look in hosts
idx = find(cellfun(@any,strfind(hostnames,name)));
if idx>0
    node = self.hosts{idx};
    return
end
% look in tasks
idx = find(cellfun(@any,strfind(tasknames,name)));
if idx>0
    node = self.tasks{idx};
    return
end
% look in entries
idx = find(cellfun(@any,strfind(entrynames,name)));
if idx>0
    node = self.entries{idx};
    return
end
% look in activities
idx = find(cellfun(@any,strfind(actnames,name)));
if idx>0
    node = self.activities{idx};
    return
end
node = [];
end
