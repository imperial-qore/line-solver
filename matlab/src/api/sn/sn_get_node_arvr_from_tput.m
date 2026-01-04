%{ @file sn_get_node_arvr_from_tput.m
 %  @brief Computes average arrival rates at nodes from station throughputs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes average arrival rates at nodes from station throughputs
 %
 % @details
 % This function calculates the average arrival rate at each node in steady-state
 % from the station throughputs, accounting for node-level routing and visits.
 %
 % @par Syntax:
 % @code
 % ANn = sn_get_node_arvr_from_tput(sn, TN, TH)
 % ANn = sn_get_node_arvr_from_tput(sn, TN, TH, AN)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>TN<td>Average throughputs at stations
 % <tr><td>TH<td>Throughput handles
 % <tr><td>AN<td>(Optional) Average arrival rates at stations; computed if missing
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>ANn<td>Average arrival rates at nodes
 % </table>
%}
function ANn=sn_get_node_arvr_from_tput(sn, TN, TH, AN)

I = sn.nnodes;
C = sn.nchains;
M = sn.nstations;
R = sn.nclasses;

if nargin<4
    AN = sn_get_arvr_from_tput(sn, TN, TH);
end

ANn = zeros(I,R);
if ~isempty(TH) && ~isempty(TN)
    for ist=1:M
        ind = sn.stationToNode(ist);
        ANn(ind,:) = AN(ist,:);
    end
    for ind=1:I
        for c = 1:C
            inchain = sn.inchain{c};
            refstat = sn.refstat(c);
            for r = inchain
                if sn.nodetype(ind) ~= NodeType.Source
                    switch sn.nodetype(ind)
                        case NodeType.Cache
                            % For cache nodes, only the requesting class arrives
                            % Hit/miss classes don't arrive - they leave
                            hitclass = sn.nodeparam{ind}.hitclass;
                            missclass = sn.nodeparam{ind}.missclass;
                            if ~(any(find(r==hitclass)) || any(find(r==missclass)))
                                % This is a requesting class (e.g., ClientClass)
                                % Arrival rate = total throughput for this class
                                ANn(ind, r) = (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                            end
                            % Hit/miss classes have 0 arrival rate at cache (they only depart)
                        otherwise
                            % For station nodes, the arrival rate is already set from AN
                            % Only use nodevisits for non-station nodes (like ClassSwitch, Sink)
                            % Note: nodeToStation is NaN for non-stations like Sink, so use ~(>0) check
                            if ~(sn.nodeToStation(ind) > 0)
                                ANn(ind, r) = (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * sum(TN(refstat,inchain));
                            end
                    end
                end
            end
        end
    end
    ANn(isnan(ANn)) = 0;
else
    ANn = [];
end
end