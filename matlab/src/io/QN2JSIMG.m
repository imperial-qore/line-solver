function outputFileName = QN2JSIMG(model, outputFileName, options)
% QN2JSIMG Writes a Network model to JMT JSIMG format
%
%   OUTPUTFILENAME = QN2JSIMG(MODEL) exports the Network MODEL to a
%   temporary JSIMG file and returns the file path.
%
%   OUTPUTFILENAME = QN2JSIMG(MODEL, OUTPUTFILENAME) exports to the
%   specified file path.
%
%   OUTPUTFILENAME = QN2JSIMG(MODEL, OUTPUTFILENAME, OPTIONS) uses the
%   specified solver options for export configuration.
%
% Parameters:
%   model - Network model to export
%   outputFileName - Optional output file path (default: temp file)
%   options - Optional SolverOptions with simulation parameters
%
% Returns:
%   outputFileName - Path to the created JSIMG file
%
% Example:
%   model = Network('example');
%   % ... define model ...
%   fname = QN2JSIMG(model);
%   jsimgView(fname);  % View in JMT
%
% See also: JMTIO, jsimgView
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    outputFileName = [];
end

if nargin < 3
    options = [];
end

% Create JMTIO instance with model and options
jmtio = JMTIO(model, options);

% Get network structure
sn = model.getStruct(true);

% Delegate to JMTIO writeJSIM
outputFileName = jmtio.writeJSIM(sn, outputFileName);
end
