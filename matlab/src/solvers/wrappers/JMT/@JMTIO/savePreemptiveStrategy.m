function [simDoc, section] = savePreemptiveStrategy(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEPREEMPTIVESTRATEGY(SIMDOC, SECTION, CURRENTNODE)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
visitsNode = simDoc.createElement('parameter');
visitsNode.setAttribute('array', 'true');
visitsNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategy');
visitsNode.setAttribute('name', 'PSStrategy');


sn = self.getStruct;
numOfClasses = sn.nclasses;
exportClasses = self.getExportableClasses();
i = sn.nodeToStation(ind);
for r=1:numOfClasses
    % Skip classes that should not be exported to JMT
    if ~exportClasses(r)
        continue;
    end

    refClassNode = simDoc.createElement('refClass');
    refClassNode.appendChild(simDoc.createTextNode(sn.classnames{r}));
    visitsNode.appendChild(refClassNode);
    
    subParameterNode = simDoc.createElement('subParameter');
    switch sn.sched(i)
        case SchedStrategy.PS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy');
            subParameterNode.setAttribute('name', 'EPSStrategy');
        case SchedStrategy.DPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategy');
            subParameterNode.setAttribute('name', 'DPSStrategy');
        case SchedStrategy.GPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategy');
            subParameterNode.setAttribute('name', 'GPSStrategy');
        case SchedStrategy.LPS
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategy');
            subParameterNode.setAttribute('name', 'EPSStrategy');
        case SchedStrategy.PSPRIO
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.EPSStrategyPriority');
            subParameterNode.setAttribute('name', 'EPSStrategyPriority');
        case SchedStrategy.DPSPRIO
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.DPSStrategyPriority');
            subParameterNode.setAttribute('name', 'DPSStrategyPriority');
        case SchedStrategy.GPSPRIO
            subParameterNode.setAttribute('classPath', 'jmt.engine.NetStrategies.PSStrategies.GPSStrategyPriority');
            subParameterNode.setAttribute('name', 'GPSStrategyPriority');
    end
    
    visitsNode.appendChild(subParameterNode);
    section.appendChild(visitsNode);
end
end

