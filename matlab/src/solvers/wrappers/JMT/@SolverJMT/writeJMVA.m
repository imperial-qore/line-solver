function outputFileName = writeJMVA(sn, outputFileName, options)
% FNAME = WRITEJMVA(SN, FNAME, OPTIONS)
%
% Write the chain-aggregated model to a JMVA XML file. The writer itself
% lives in @JMTIO/writeJMVA.m, next to the other JMT serializers, and is
% delegated to here rather than duplicated: the two copies had already
% drifted, with this one supporting load-dependent open models that the
% other rejected.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outputFileName = JMTIO.writeJMVA(sn, outputFileName, options);
end
