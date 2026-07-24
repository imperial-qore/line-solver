function plot(self, options)
% PLOT Display network diagram using TikZ/LaTeX
%
% PLOT() displays the network as a TikZ diagram in the system's
% default PDF viewer. This is an alias for tikzView().
%
% PLOT(OPTIONS) displays with custom TikZOptions configuration.
%
% Requires pdflatex to be installed.
%
% Example:
%   model.plot()
%   opts = jline.io.tikz.TikZOptions();
%   opts.setNodeSpacing(4.0);
%   model.plot(opts)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    self.tikzView();
else
    self.tikzView(options);
end
end
