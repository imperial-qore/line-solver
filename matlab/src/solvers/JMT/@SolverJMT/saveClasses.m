function [simElem, simDoc] = saveClasses(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVECLASSES(SIMELEM, SIMDOC)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sn = self.getStruct;

numOfClasses = sn.nclasses;
for r=1:numOfClasses
    userClass = simDoc.createElement('userClass');
    userClass.setAttribute('name', sn.classnames(r));
    if isinf(sn.njobs(r))
        userClass.setAttribute('type', 'open');
    else
        userClass.setAttribute('type', 'closed');
    end
    userClass.setAttribute('priority', int2str(sn.classprio(r)));
    refStatIndex = sn.refstat(r);
    refNodeIndex = sn.stationToNode(sn.refstat(r));
    refStatName = sn.nodenames{refNodeIndex};
    if ~isempty(sn.proc{refStatIndex}{r})
        if isfinite(sn.njobs(r)) % if closed
            userClass.setAttribute('customers', int2str(sn.njobs(r)));
            userClass.setAttribute('referenceSource', refStatName);
        elseif isnan(sn.proc{refStatIndex}{r}{1}) % open disabled in source
            userClass.setAttribute('referenceSource', 'ClassSwitch');
        else % if other open
            userClass.setAttribute('referenceSource', sn.nodenames{sn.stationToNode(sn.refstat(r))});
        end
    else
        userClass.setAttribute('referenceSource', sn.nodenames{sn.stationToNode(sn.refstat(r))});
    end
    simElem.appendChild(userClass);
end

end
