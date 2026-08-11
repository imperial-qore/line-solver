function entries = citations(self)
% ENTRIES = CITATIONS()
%
% Bibliographic references for the solver-selection policy this run used, as
% a struct array with fields .key, .ref and .covers, printed as a list when
% called without an output argument.
%
% SolverAUTO performs no numerical work of its own, so it cites only how the
% solver was chosen. The algorithms that produced the numbers belong to the
% delegate: ask it with getSolver().citations() when the run has one. The
% default heuristic policy is hand-written and has no reference; the learned
% policy ('method','tree') does.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tokens = {};
if strcmp(self.selectionMode, 'tree')
    tokens{end+1} = 'auto.tree';
    tokens{end+1} = 'cart';
end

entries = line_citations(tokens);

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
