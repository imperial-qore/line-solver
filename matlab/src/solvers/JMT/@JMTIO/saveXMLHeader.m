function [simElem,simDoc] = saveXMLHeader(self, logPath)
% [SIMELEM,SIMDOC] = SAVEXMLHEADER(LOGPATH)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
xmlnsXsi = 'http://www.w3.org/2001/XMLSchema-instance';
fname = [getFileName(self), ['.', 'jsimg']];
simDoc = com.mathworks.xml.XMLUtils.createDocument('sim');
simElem = simDoc.getDocumentElement;
simElem.setAttribute('xmlns:xsi', xmlnsXsi);
simElem.setAttribute('name', fname);
%simElem.setAttribute('timestamp', '"Tue Jan 1 00:00:01 GMT+00:00 2000"');
simElem.setAttribute('xsi:noNamespaceSchemaLocation', 'SIMmodeldefinition.xsd');
simElem.setAttribute('disableStatisticStop', 'true');
simElem.setAttribute('logDecimalSeparator', '.');
simElem.setAttribute('logDelimiter', ';');
simElem.setAttribute('logPath', logPath);
simElem.setAttribute('logReplaceMode', '0');
% Set maxSamples (default to 10000 if not set)
if isempty(self.maxSamples) || isnan(self.maxSamples)
    simElem.setAttribute('maxSamples', '10000');
else
    simElem.setAttribute('maxSamples', int2str(self.maxSamples));
end
simElem.setAttribute('maxEvents', int2str(self.maxEvents));
if ~isinf(self.maxSimulatedTime)
    simElem.setAttribute('maxSimulated', num2str(self.maxSimulatedTime,'%.3f'));
end
simElem.setAttribute('polling', '1.0');
simElem.setAttribute('seed', int2str(self.options.seed));
end
