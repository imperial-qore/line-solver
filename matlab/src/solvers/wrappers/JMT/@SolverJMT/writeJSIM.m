function outputFileName = writeJSIM(self, sn, outputFileName)
% FNAME = WRITEJSIM(SN, FNAME)
%
% Delegates XML generation to JMTXMLParser.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(sn)
    sn = self.getStruct();
end

if nargin < 3 || isempty(outputFileName)
    outputFileName = getJSIMTempPath(self);
end

% Pass handles to JMTXMLParser
self.xmlParser.handles = self.handles;

% Delegate to JMTXMLParser
outputFileName = self.xmlParser.writeJSIM(sn, outputFileName);
end
