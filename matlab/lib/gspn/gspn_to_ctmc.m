function [State, Adj, StateTable] = gspn_to_ctmc(F, B, INH, IM, P, M0)
% GSPN_TO_CTMC  Generate the CTMC from a gspn
% F - Forward Matrix
% B - Backward Matrix
% INH - Matrix that represents inhibitor arcs
% IM - Boolean matrix that distinguish between immediate and timed
% transitions
% P - Matrix that represents priorities between immediate transitions.
% W - Matrix that represents weights between immediate transitions.
% D - Matrix that holds the distribution of different modes of transition.

% nt - number of transitions
[num_m, num_t] = size(F);
State = M0;
StateTable = HashTable(97);
StateTable.insert(M0);
Adj = AdjacencyList(97);

C = cellfun(@minus,B,F,'Un',0);
i=0;

while i < size(State,2)
    i=i+1;
    enabled_t = [];
    enabled_m = [];
    for t=1:num_t
        for m=1:num_m
            fi = cell2mat(F(m,t));
            if ~isempty(fi)
                in = cell2mat(INH(m,t));
                enabled = all(State(:,i) >= fi & ~any(in>0 & in <= State(:,i)));
                if enabled
                    if isempty(enabled_m)
                        enabled_m(end+1) = m;
                        enabled_t(end+1) = t;
                    else
                        em = enabled_m(end);
                        et = enabled_t(end);
                        isImmed_1 = cell2mat(IM(m,t));
                        isImmed_2 = cell2mat(IM(em,et));
                        if isImmed_1 && isImmed_2
                            % Both mode is immediate
                            priority_1 = cell2mat(P(m,t));
                            priority_2 = cell2mat(P(em,et));
                            if priority_1 > priority_2
                                enabled_m = m;
                                enabled_t = t;
                            elseif priority_1 == priority_2
                                enabled_m(end+1) = m;
                                enabled_t(end+1) = t;
                            end
                        elseif isImmed_1 && ~isImmed_2
                            enabled_m = m;
                            enabled_t = t;
                        elseif ~isImmed_1 && ~isImmed_2
                            enabled_m(end+1) = m;
                            enabled_t(end+1) = t;
                        end
                    end
                end
            end
        end
    end
    for j=1:length(enabled_m)
        cur_c = cell2mat(C(enabled_m(j),enabled_t(j)));
        new_state = State(:,i)+cur_c;
        [exists,sequence_number] = StateTable.search(new_state);
        if exists
%             sprintf('Duplicate State Found')
            Adj.insert(i, sequence_number, [enabled_m(j),enabled_t(j)]);
        else
            State=[State,new_state];
            sequence_number = StateTable.insert(new_state);
            Adj.insert(i, sequence_number, [enabled_m(j),enabled_t(j)]);
        end
    end
end