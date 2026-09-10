%{ @file sn_pn_avg_rates.m
 %  @brief Converts Place throughputs from firing events to tokens and derives
 %         the matching arrival and response times
 %
 %  @author LINE Development Team
%}

%{
 % @brief Place throughput, arrival rate and response time in tokens
 %
 % @details
 % A Place is a station and a token is the job it holds, so a firing that
 % consumes two tokens is two departures, not one. The CTMC and SSA analyzers
 % count firing events instead, which for unit arc multiplicities is the same
 % number and for weighted arcs is not: the reported throughput is then not a
 % token rate, and QLen over it is not a sojourn time. On a net where a Place
 % is drained by an arc of weight 2 and another of weight 3, the event count
 % gave RespT 1.6243 where Little's law on tokens gives 0.7119, which is the
 % value SolverJMT measures.
 %
 % This function rescales the Place rows to tokens:
 %
 %   TN(p,k)  tokens consumed from the Place per unit time
 %   AN(p,k)  tokens produced into the Place per unit time
 %   RN(p,k)  QN(p,k) / TN(p,k), Little's law over the Place
 %
 % Rows that do not belong to a Place are returned untouched, so a mixed
 % Queue/Place model keeps its queueing metrics. When the firing rates cannot
 % be recovered from the throughputs the inputs are returned unchanged rather
 % than replaced by a guess.
 %
 % @par Syntax:
 % @code
 % [TN, AN, RN] = sn_pn_avg_rates(sn, QN, TN, AN, RN)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>QN<td>Average queue lengths, i.e. mean token counts at the Places
 % <tr><td>TN<td>Average throughputs at stations, counting firing events
 % <tr><td>AN<td>Average arrival rates at stations, as computed by the caller
 % <tr><td>RN<td>Average response times at stations, as computed by the caller
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>TN<td>Throughputs, with the Place rows counting tokens
 % <tr><td>AN<td>Arrival rates, with the Place rows counting tokens
 % <tr><td>RN<td>Response times, with the Place rows following Little's law
 % </table>
%}
function [TN, AN, RN] = sn_pn_avg_rates(sn, QN, TN, AN, RN)

if isempty(TN) || ~any(sn.nodetype == NodeType.Place)
    return
end

% The analyzers hand over event counts, hence the false.
[x, consumed, produced, placeNodes] = sn_pn_firing_rates(sn, TN, false);
if isempty(x)
    return
end

R = sn.nclasses;
for pp = 1:length(placeNodes)
    ist = sn.nodeToStation(placeNodes(pp));
    if ~(ist > 0)
        continue
    end
    for k = 1:R
        TN(ist, k) = consumed(:, pp, k)' * x;
        if ~isempty(AN)
            AN(ist, k) = produced(:, pp, k)' * x;
        end
        if ~isempty(RN)
            if TN(ist, k) > 0
                RN(ist, k) = QN(ist, k) / TN(ist, k);
            else
                RN(ist, k) = 0;
            end
        end
    end
end
end
