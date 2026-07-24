function outputFileName = QN2JSIMG(self, sn, outputFileName)
% FNAME = QN2JSIMG(SN, FNAME)
%
% Writes queueing network model to JMT JSIMG format.
% Delegates to the standalone QN2JSIMG function in the io package.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3
    outputFileName = [];
end

% Delegate to io package function
outputFileName = QN2JSIMG(self.model, outputFileName, self.options);
end
