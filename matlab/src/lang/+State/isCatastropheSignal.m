function bool = isCatastropheSignal(sn, class)
% BOOL = ISCATASTROPHESIGNAL(SN, CLASS)
%
% True if signal class CLASS empties a station on arrival. The signaltype is
% consulted alongside the iscatastrophe flag so that the two encodings of a
% catastrophe cannot disagree (SolverMAM applies the same test).

bool = false;
if isfield(sn, 'iscatastrophe') && ~isempty(sn.iscatastrophe) && numel(sn.iscatastrophe) >= class ...
        && sn.iscatastrophe(class)
    bool = true;
    return
end
if isfield(sn, 'signaltype') && ~isempty(sn.signaltype) && numel(sn.signaltype) >= class
    st = sn.signaltype{class};
    if ~isempty(st) && ~(isnumeric(st) && any(isnan(st))) && st == SignalType.CATASTROPHE
        bool = true;
    end
end
end
