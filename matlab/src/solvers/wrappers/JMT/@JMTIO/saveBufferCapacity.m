function [simDoc, section] = saveBufferCapacity(self, simDoc, section, ind)
% [SIMDOC, SECTION] = SAVEBUFFERCAPACITY(SIMDOC, SECTION, NODEIDX)
%
% LINE uses Kendall notation where cap = K = total system capacity
% (queue + in-service jobs). JMT's "size" parameter also represents
% total capacity K.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.getStruct;
sizeNode = simDoc.createElement('parameter');
sizeNode.setAttribute('classPath', 'java.lang.Integer');
sizeNode.setAttribute('name', 'size');
valueNode = simDoc.createElement('value');
ist = sn.nodeToStation(ind);
if ~sn.isstation(ind) || isinf(sn.cap(ist))
    valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
else
    if sn.cap(ist) >= jmtReachablePopulation(sn, ist)
        % A capacity the population cannot reach is unbounded, and the test is
        % >= and not ==: refreshCapacity DERIVES sn.cap for a station the user
        % never capped, as sum over the classes served there of the chain
        % population, so a multi-class station gets (#classes) x N -- 8 on a
        % two-class model of 4 jobs. Under == only the single-class case
        % matched, and every multi-class one fell through to
        % jmtStationCapRefusal and was refused as a "binding" buffer nobody
        % declared.
        %
        % THE POPULATION COMPARED AGAINST IS THE ONE THAT CAN REACH IST, not
        % sum(sn.njobs). A class that never visits this station cannot fill it,
        % so counting its jobs makes a capacity that is exactly the reachable
        % population look like a buffer. That is what a SELF-LOOPING CLASS does:
        % on sanity_CQN_rm_fcfs_1class_2slcatdelay the SLC's 2 jobs stay at the
        % Delay, Queue1 derives cap = 1 for the one class that visits it, and
        % 1 >= 3 is false -- so SolverJMT refused, with a message naming a
        % finite capacity the model never declares, a model it had always
        % exported. Same on sanity_CQN_2q_psfcfs_1class_1slcateachqueue.
        valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
    else
        nservers = sn.nservers(ist);
        if isinf(nservers)
            % Infinite servers (delay node) - no buffer needed
            valueNode.appendChild(simDoc.createTextNode(int2str(-1)));
        else
            % The same predicate jmtMethodRefusal asks, so the gate that
            % decides whether to OFFER jmt.jsim and this writer cannot answer
            % differently. It was a local function here, which is exactly why
            % the gate could not see it and offered a pair that then raised.
            capReason = jmtStationCapRefusal(sn, ist);
            if ~isempty(capReason)
                line_error(mfilename, capReason);
            end
            % Send LINE's total capacity K directly to JMT
            valueNode.appendChild(simDoc.createTextNode(int2str(sn.cap(ist))));
        end
    end
end

sizeNode.appendChild(valueNode);
section.appendChild(sizeNode);
end
