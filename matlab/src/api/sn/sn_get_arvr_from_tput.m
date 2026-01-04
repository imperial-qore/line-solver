%{ @file sn_get_arvr_from_tput.m
 %  @brief Computes average arrival rates at stations from throughputs
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes average arrival rates at stations from throughputs
 %
 % @details
 % This function calculates the average arrival rate at each station
 % in steady-state from the station throughputs and routing matrix.
 %
 % @par Syntax:
 % @code
 % AN = sn_get_arvr_from_tput(sn, TN, TH)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>TN<td>Average throughputs at stations
 % <tr><td>TH<td>Throughput handles (optional)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>AN<td>Average arrival rates at stations
 % </table>
%}
function AN=sn_get_arvr_from_tput(sn, TN, TH)

M = sn.nstations;
R = sn.nclasses;
if ~isempty(TH) && ~isempty(TN)
    AN = zeros(M,R);

    % Build mapping from stateful nodes to their position in rt matrix
    % rt is indexed by stateful nodes, not stations
    statefulNodes = find(sn.isstateful);
    nStateful = length(statefulNodes);

    % Build throughput vector for all stateful nodes (stations have TN, others need computation)
    TN_stateful = zeros(nStateful, R);
    for sf = 1:nStateful
        ind = statefulNodes(sf);
        ist = sn.nodeToStation(ind);
        if ist > 0
            % This stateful node is a station - use station throughput
            TN_stateful(sf, :) = TN(ist, :);
        else
            % This stateful node is not a station (e.g., Cache)
            % Compute throughput from routing and reference station throughput
            if sn.nodetype(ind) == NodeType.Cache
                % For Cache nodes, compute hit/miss class throughputs
                % from the reference station throughput and hit/miss probabilities
                hitclass = sn.nodeparam{ind}.hitclass;
                missclass = sn.nodeparam{ind}.missclass;

                % Get actual hit/miss probabilities if available
                if isfield(sn.nodeparam{ind}, 'actualhitprob') && ~isempty(sn.nodeparam{ind}.actualhitprob)
                    actualHitProb = sn.nodeparam{ind}.actualhitprob;
                    actualMissProb = sn.nodeparam{ind}.actualmissprob;
                else
                    % Actual probabilities not yet computed - skip this cache
                    % Arrival rates will be computed later when probabilities are available
                    continue;
                end

                % Find the chain and reference station for this cache
                for c = 1:sn.nchains
                    inchain = sn.inchain{c};
                    refstat = sn.refstat(c);
                    totalTput = sum(TN(refstat, inchain));

                    % Set throughput for hit/miss classes
                    for origClass = 1:length(hitclass)
                        if hitclass(origClass) > 0 && hitclass(origClass) <= R
                            TN_stateful(sf, hitclass(origClass)) = totalTput * actualHitProb(origClass);
                        end
                        if missclass(origClass) > 0 && missclass(origClass) <= R
                            TN_stateful(sf, missclass(origClass)) = totalTput * actualMissProb(origClass);
                        end
                    end
                end
            end
        end
    end

    % Compute arrival rates using stateful node throughputs and rt matrix
    for ist=1:M
        ind_ist = sn.stationToNode(ist);
        if sn.nodetype(ind_ist)==NodeType.Source
            AN(ist,:) = 0;
        else
            % Find the position of this station in the stateful node list
            sf_ist = find(statefulNodes == ind_ist);
            if isempty(sf_ist)
                continue;
            end

            for sf_jst = 1:nStateful
                for k=1:R
                    for r=1:R
                        AN(ist,k) = AN(ist,k) + TN_stateful(sf_jst,r)*sn.rt((sf_jst-1)*R+r, (sf_ist-1)*R+k);
                    end
                end
            end
        end
    end
else
    AN = [];
end

if any(sn.fj(:))
    ANn = sn_get_node_arvr_from_tput(sn, TN, TH, AN);
    for ist=1:M
        AN(ist,:) = ANn(sn.stationToNode(ist),:);
    end
end

end