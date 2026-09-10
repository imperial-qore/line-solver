function entries = citations(self)
% ENTRIES = CITATIONS()
%
% Bibliographic references for the solver-selection policy this run used, as
% a struct array with fields .key, .ref and .covers, printed as a list when
% called without an output argument.
%
% SolverAUTO performs no numerical work of its own, so it cites only how the
% solver was chosen. The algorithms that produced the numbers belong to the
% delegate: ask it with getSolver().citations() when the run has one. Every
% selection policy is hand-written and has no reference, so the list is empty.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

methodNames = {};

entries = line_citations(methodNames);

if nargout == 0
    if isempty(entries)
        line_printf('No selection-policy references recorded for this run.\n');
    else
        for i = 1:numel(entries)
            line_printf('%s\n    covers: %s\n', entries(i).ref, entries(i).covers);
        end
    end
    clear entries
end
end
