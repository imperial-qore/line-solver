function view(self, options)
% VIEW(OPTIONS) Alias for jsimgView
%
% Opens the model in JMT's graphical editor (JSIMgraph).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = Solver.defaultOptions();
end
self.jsimgView(options);
end
