function [simDoc, section] = savePreemptiveStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEPREEMPTIVESTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
visitsNode = simDoc.createElement('parameter');
visitsNode.setAttribute('array', 'true');
visitsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategy');
visitsNode.setAttribute('name', 'PSStrategy');


sn = self.getStruct;
numOfClasses = sn.nclasses;
i = sn.nodeToStation(ind);
for r=1:numOfClasses  
    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    visitsNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    switch sn.schedid(i)
        case SchedStrategy.ID_PS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy');
            subParameterNode.setAttribute('name', 'EPSStrategy');
        case SchedStrategy.ID_DPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategy');
            subParameterNode.setAttribute('name', 'DPSStrategy');
        case SchedStrategy.ID_GPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategy');
            subParameterNode.setAttribute('name', 'GPSStrategy');
    end
    
    visitsNode.appendChild(subParameterNode);
    section.appendChild(visitsNode);
end
end

