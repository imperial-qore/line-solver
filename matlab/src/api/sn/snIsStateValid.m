function [isvalid] = snIsStateValid(sn)
% [ISVALID] = ISSTATEVALID()
nir = [];
sir = [];
for ist=1:sn.nstations
    isf = sn.stationToStateful(ist);
    if size(sn.state{isf},1)>1
        line_warning(mfilename,sprintf('isStateValid will ignore some states of station %d, define a unique initial state to address this problem.\n',ist));
        sn.state{isf} = sn.state{isf}(1,:);
    end
    [~, nir(ist,:), sir(ist,:), ~] = State.toMarginal(sn, sn.stationToNode(ist), sn.state{isf});
end
isvalid = State.isValid(sn, nir, sir);
end