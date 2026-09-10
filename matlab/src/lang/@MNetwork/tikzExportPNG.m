function tikzExportPNG(self, filePath, dpi)
% TIKZEXPORTPNG Export network diagram to PNG image
%
% TIKZEXPORTPNG(FILEPATH) exports the network diagram as a PNG image
% at 150 DPI resolution.
%
% TIKZEXPORTPNG(FILEPATH, DPI) exports at the specified DPI resolution.
%
% Requires pdflatex and pdftoppm (from poppler-utils) to be installed.
%
% Example:
%   model.tikzExportPNG('network.png')
%   model.tikzExportPNG('network_hires.png', 300)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Convert MATLAB model to Java Network
jnetwork = JLINE.line_to_jline(self);

if nargin < 3
    jnetwork.tikzExportPNG(filePath);
else
    jnetwork.tikzExportPNG(filePath, int32(dpi));
end
end
