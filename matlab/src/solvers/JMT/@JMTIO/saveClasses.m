function [simElem, simDoc] = saveClasses(self, simElem, simDoc)
% [SIMELEM, SIMDOC] = SAVECLASSES(SIMELEM, SIMDOC)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;

% Get exportable classes (handles cache classes and class-switching)
exportClasses = self.getExportableClasses();

% JMT uses higher priority value = higher priority, while LINE uses lower value = higher priority
% We need to invert priorities when writing to JMT
maxPrio = max(sn.classprio);

numOfClasses = sn.nclasses;
for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    userClass = simDoc.createElement('userClass');
    userClass.setAttribute('name', sn.classnames(r));
    if isinf(sn.njobs(r))
        userClass.setAttribute('type', 'open');
    else
        userClass.setAttribute('type', 'closed');
    end
    % Set soft deadline from class deadline (0.0 if no deadline)
    deadline = sn.classdeadline(r);
    if isfinite(deadline)
        userClass.setAttribute('softDeadline', num2str(deadline));
    else
        userClass.setAttribute('softDeadline', '0.0');
    end
    % Invert priority: LINE uses lower=higher, JMT uses higher=higher
    jmtPrio = maxPrio - sn.classprio(r);
    userClass.setAttribute('priority', int2str(jmtPrio));
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
