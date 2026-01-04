function self = relink(self, P)
% SELF = RELINK(P)
% Call LINK again with a new P

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

self.resetNetwork();
self = link(self, P);
end
