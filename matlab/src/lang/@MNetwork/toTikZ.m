function tikzCode = toTikZ(self, options)
% TOTIKZ Generate TikZ/LaTeX code for the network diagram
%
% CODE = TOTIKZ() returns a complete LaTeX document containing
% the TikZ code that can be compiled with pdflatex.
%
% CODE = TOTIKZ(OPTIONS) generates with custom TikZOptions configuration.
%
% Example:
%   tikzCode = model.toTikZ();
%   fid = fopen('network.tex', 'w');
%   fprintf(fid, '%s', tikzCode);
%   fclose(fid);
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Convert MATLAB model to Java Network
jnetwork = JLINE.line_to_jline(self);

if nargin < 2
    tikzCode = char(jnetwork.toTikZ());
else
    tikzCode = char(jnetwork.toTikZ(options));
end
end
