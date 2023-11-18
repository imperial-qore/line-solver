function [ni, nir, sir, kir] = snInitToMarginal(sn)
% [NI, NIR, SIR, KIR] = INITTOMARGINAL()

ni = {}; nir = {}; sir = {}; kir = {};
for ist=1:sn.nstations
    if ~isempty(sn.state{ist})
        [ni{ist,1}, nir{ist,1}, sir{ist,1}, kir{ist,1}] = State.toMarginal(sn,sn.stationToNode(ist),state);
    end
end
end