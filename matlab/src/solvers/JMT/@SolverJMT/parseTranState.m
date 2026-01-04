function [state, evtype, evclass, evjob] = parseTranState(fileArv, fileDep, nodePreload)
% [STATE, EVTYPE, EVCLASS, EVJOB] = PARSETRANSTATE(FILEARV, FILEDEP, NODEPRELOAD)
%
% Delegates to JMTIO.parseTranState for backward compatibility.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[state, evtype, evclass, evjob] = JMTIO.parseTranState(fileArv, fileDep, nodePreload);
end
