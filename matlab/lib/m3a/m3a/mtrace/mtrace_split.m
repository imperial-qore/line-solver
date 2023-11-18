function [TL] = mtrace_split(T, L)
% Given a multi-class trace with inter-arrivals T and labels L,
% creates the separate per-class traces. In each per-class trace
% the inter-arrivals are inter-arrivals between events of the same class.
% The result is a cell array and the c-th element of the cell array
% contains the vector of inter-arrival times for class c.

labels = unique(L);
C = length(labels);

TL = cell(1,C);

TCUM = cumsum(T);
for c = 1:C
    TL{c} = diff([0; TCUM(L == labels(c))]);
end

end