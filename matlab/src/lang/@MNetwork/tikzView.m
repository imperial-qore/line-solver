function tikzView(self, options)
% TIKZVIEW Display network diagram using TikZ/LaTeX
%
% TIKZVIEW() displays the network as a TikZ diagram in the system's
% default PDF viewer. Requires pdflatex to be installed.
%
% TIKZVIEW(OPTIONS) displays with custom TikZOptions configuration.
%
% Example:
%   model.tikzView()
%   opts = jline.io.tikz.TikZOptions();
%   opts.setShowRoutingProb(false);
%   model.tikzView(opts)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Convert MATLAB model to Java Network
jnetwork = JLINE.line_to_jline(self);

if nargin < 2
    jnetwork.tikzView();
else
    jnetwork.tikzView(options);
end
end
