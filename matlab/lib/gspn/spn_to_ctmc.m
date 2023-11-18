function [State, Adj] = spn_to_ctmc(F, B, INH, M0)

% nt - number of transitions
[~, num_t] = size(F);
State = M0;
Adj = [0];
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
        new_state_rep = repmat(new_state,1,size(State,2));
        duplicates =all(new_state_rep == State);
		duplicate_found =find(duplicates);
        if size(duplicate_found,2)>1
            sprintf('Duplicate State Found')
        end
        if any(duplicates)
            Adj(i,duplicate_found)=enabled_transitions(t);
        else
			State=[State,new_state];
			Adj(length(Adj)+1,size(Adj,2)+1)=0;
			Adj(i,size(Adj,2))=enabled_transitions(t);
        end
    end
end