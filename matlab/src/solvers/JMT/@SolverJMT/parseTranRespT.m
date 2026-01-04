function [classResT, jobResT, jobResTArvTS, classResTJobID] = parseTranRespT(fileArv, fileDep)
% [CLASSREST, JOBREST, JOBRESTARTVTS, CLASSRESTJOBID] = PARSETRANRESPT(FILEARV, FILEDEP)
%
% Delegates to JMTIO.parseTranRespT for backward compatibility.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[classResT, jobResT, jobResTArvTS, classResTJobID] = JMTIO.parseTranRespT(fileArv, fileDep);
end
