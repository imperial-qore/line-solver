function [State, Adj, StateTable] = spn_to_ctmc_hashtable(F, B, INH, M0)

% nt - number of transitions
[~, num_t] = size(F);
State = M0;
StateTable = HashTable(97);
StateTable.insert(M0);
Adj = AdjacencyList(97);

C = B-F;
i=0;

while i < size(State,2)
    i=i+1;
    
    for t=1:num_t
        x(t) = all(State(:,i) >= F(:,t) & ~any(INH(:,t)>0 & INH(:,t) <= State(:,i)));
    end
    enabled_transitions = find(x);
    
    for t=1:size(enabled_transitions,2)
        new_state = State(:,i)+C(:,enabled_transitions(t));
        [exists,sequence_number] = StateTable.search(new_state);
        if exists
            sprintf('Duplicate State Found')
            Adj.insert(i, sequence_number, enabled_transitions(t));
        else
            State=[State,new_state];
            sequence_number = StateTable.insert(new_state);
            Adj.insert(i, sequence_number, enabled_transitions(t));
        end
    end
end