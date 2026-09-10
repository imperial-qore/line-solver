function logData = parseLogs(model, isNodeLogged, metric)
% LOGDATA = PARSELOGS(MODEL,ISNODELOGGED, METRIC)
%
% Delegates to JMTIO.parseLogs for backward compatibility.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

logData = JMTIO.parseLogs(model, isNodeLogged, metric);
end
