function [simDoc, section] = saveClassSwitchStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVECLASSSWITCHSTRATEGY(SIMDOC, SECTION, NODEIDX)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

paramNode = simDoc.createElement('parameter');
paramNode.setAttribute('array', 'true');
paramNode.setAttribute('classPath', 'java.lang.Object');
paramNode.setAttribute('name', 'matrix');

sn = self.getStruct;
K = sn.nclasses;
exportClasses = self.getExportableClasses();
i = ind;
jset = find(sn.connmatrix(ind,:));
for r=1:K
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    paramNode.appendChild(refClassNode);


    subParNodeRow = simDoc.createElement('subParameter');
    subParNodeRow.setAttribute('array', 'true');
    subParNodeRow.setAttribute('classPath', 'java.lang.Float');
    subParNodeRow.setAttribute('name', 'row');
    for s=1:K
        % Skip classes that should not be exported to JMT
        if ~exportClasses(s)
            continue;
        end

        refClassNode = simDoc.createElement('refClass');
        refClassNode.appendChild(simDoc.createTextNode(sn.classnames{s}));
        subParNodeRow.appendChild(refClassNode);
        
        subParNodeCell = simDoc.createElement('subParameter');
        subParNodeCell.setAttribute('classPath', 'java.lang.Float');
        subParNodeCell.setAttribute('name', 'cell');
        valNode = simDoc.createElement('value');
        valNode.appendChild(simDoc.createTextNode(sprintf('%12.12f', sum(sn.rtnodes((i-1)*K+r, (jset-1)*K+s)))));
        subParNodeCell.appendChild(valNode);
        subParNodeRow.appendChild(subParNodeCell);
        
    end
    paramNode.appendChild(subParNodeRow);
    
end
section.appendChild(paramNode);

end
