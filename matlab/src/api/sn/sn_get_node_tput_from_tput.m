%{ @file sn_get_node_tput_from_tput.m
 %  @brief Computes average throughputs at nodes from station throughputs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes average throughputs at nodes from station throughputs
 %
 % @details
 % This function calculates the average throughput at each node in steady-state
 % from the station throughputs and node-level routing matrix.
 %
 % @par Syntax:
 % @code
 % TNn = sn_get_node_tput_from_tput(sn, TN, TH)
 % TNn = sn_get_node_tput_from_tput(sn, TN, TH, ANn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>TN<td>Average throughputs at stations
 % <tr><td>TH<td>Throughput handles
 % <tr><td>ANn<td>(Optional) Average arrival rates at nodes; computed if missing
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>TNn<td>Average throughputs at nodes
 % </table>
%}
function TNn=sn_get_node_tput_from_tput(sn, TN, TH, ANn)

I = sn.nnodes;
C = sn.nchains;
R = sn.nclasses;
if nargin<4
    ANn = sn_get_node_arvr_from_tput(sn, TN, TH);
end

TNn = zeros(I,R);
if ~isempty(TH) && ~isempty(TN)
    for ind=1:I
        for c = 1:C
            inchain = sn.inchain{c};
            refstat = sn.refstat(c);
            for r = inchain
                if sn.nodetype(ind) ~= NodeType.Source
                    switch sn.nodetype(ind)
                        case NodeType.Cache
                            % For cache nodes, use actual hit/miss ratios if available
                            % instead of nodevisits which don't account for cache behavior
                            hitclass = sn.nodeparam{ind}.hitclass;
                            missclass = sn.nodeparam{ind}.missclass;
                            totalTput = sum(TN(refstat,inchain));

                            % Check if actual hit/miss probabilities are available in nodeparam
                            if isfield(sn.nodeparam{ind}, 'actualhitprob') && ~isempty(sn.nodeparam{ind}.actualhitprob)
                                % Use actual hit/miss probabilities from solver result
                                actualHitProb = sn.nodeparam{ind}.actualhitprob;
                                actualMissProb = sn.nodeparam{ind}.actualmissprob;

                                % Find which original class this hit/miss class belongs to
                                for origClass = 1:length(hitclass)
                                    if hitclass(origClass) == r
                                        % This is a hit class - use hit probability
                                        TNn(ind, r) = totalTput * actualHitProb(origClass);
                                    elseif missclass(origClass) == r
                                        % This is a miss class - use miss probability
                                        TNn(ind, r) = totalTput * actualMissProb(origClass);
                                    end
                                end
                            else
                                % Fallback to nodevisits-based calculation
                                if any(find(r==hitclass)) || any(find(r==missclass))
                                    TNn(ind, r) = (sn.nodevisits{c}(ind,r) / sum(sn.visits{c}(sn.stationToStateful(refstat),inchain))) * totalTput;
                                end
                            end
                    end
                end
            end
        end
    end

    % First, copy station throughputs directly to station nodes
    M = sn.nstations;
    for ist=1:M
        ind = sn.stationToNode(ist);
        TNn(ind,:) = TN(ist,:);
    end

    for ind=1:I
        for c = 1:C
            inchain = sn.inchain{c};
            for r = inchain
                anystateful = find(sn.visits{c}(:,r));
                if ~isempty(anystateful)
                    if sn.nodetype(ind) ~= NodeType.Sink && sn.nodetype(ind) ~= NodeType.Join
                        for s = inchain
                            for jnd=1:I
                                switch sn.nodetype(ind)
                                    case NodeType.Source
                                        ist = sn.nodeToStation(ind);
                                        TNn(ind, s) = TN(ist,s);
                                    case NodeType.Cache
                                        if ind~=jnd
                                            TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                                        end
                                    otherwise
                                        % For station nodes, throughput is already set from TN
                                        % Only compute for non-station nodes (like ClassSwitch, Router)
                                        % Note: nodeToStation returns NaN for non-station nodes, so use ~(>0) check
                                        if ~(sn.nodeToStation(ind) > 0)
                                            TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                                        end
                                end
                            end
                        end
                    elseif sn.nodetype(ind) == NodeType.Join
                        for s = inchain
                            for jnd=1:I
                                if sn.nodetype(ind) ~= NodeType.Source
                                    TNn(ind, s) = TNn(ind, s) + ANn(ind, r) * sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                                else
                                    ist = sn.nodeToStation(ind);
                                    TNn(ind, s) = TN(ist,s);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    TNn(isnan(TNn)) = 0;
else
    TNn = [];
end

end