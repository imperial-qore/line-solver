function [State, Adj, StateTable] = spn_to_ctmc_parallel(F, B, INH, M0)

% nt - number of transitions
[~, num_t] = size(F);
% workers = parpool('local');
% N = workers.NumWorkers;

spmd(8)
    N = 8;
    e = 1;
    initial_node = h0(M0, N);
    if initial_node == labindex
        State = M0;
        StateTable = HashTable(97);
        StateTable.insert(M0);
    else
        State = [];
        StateTable = HashTable(97);
    end
    
    Adj = AdjacencyList(97);
    C = B-F;
    i=0;
    
    while e < 500
        e = e + 1;
        if i < size(State,2)
            i=i+1;
            for t=1:num_t
                x(t) = all(State(:,i) >= F(:,t) & ~any(INH(:,t)>0 & INH(:,t) <= State(:,i)));
            end
            enabled_transitions = find(x);
            for t=1:size(enabled_transitions,2)
                new_state = State(:,i)+C(:,enabled_transitions(t));
                fprintf('new_state [')
                for c=1:length(new_state)
                    fprintf('%d ',new_state(c))
                end
                fprintf(']\n');
                assigned_node = h0(new_state, N);
                fprintf('To %d\n',assigned_node)
                if assigned_node == labindex
                    [exists,sequence_number] = StateTable.search(new_state);
                    if exists
                        sprintf('Duplicate State Found')
                        Adj.insert(i, [labindex, sequence_number], enabled_transitions(t));
%                         sprintf('Endssss')
                    else
                        State=[State,new_state];
                        sequence_number = StateTable.insert(new_state);
                        Adj.insert(i, [labindex, sequence_number], enabled_transitions(t));
                    end
                else
                    fprintf('Sent Data %s\n',class(new_state))
                    labSend({[labindex, i],new_state, enabled_transitions(t)}, assigned_node, 1); % Tag == 1 Send State
                end
            end
        end
        while labProbe('any', 2) % Tag == 2 IDs
%             sprintf('Receiving IdInfo')
            [idInfo, src] = labReceive('any',2);
            parentSeq = cell2mat(idInfo(1));
            childSeq = cell2mat(idInfo(2));
            enabled_transition = cell2mat(idInfo(3));
            Adj.insert(parentSeq, [src, childSeq], enabled_transition);
        end
        while labProbe('any', 1) % Tag == 2 State
%             sprintf('Receiving StateInfo')
            [stateInfo, src] = labReceive('any',1); % StateInfo contains id of parents and state descriptor
            id_parent = cell2mat(stateInfo(1)); % [nodeIdx, sequence_number]
            new_state = cell2mat(stateInfo(2));
            enabled_transition = cell2mat(stateInfo(3));
            [exists,sequence_number] = StateTable.search(new_state);
            if exists
                sprintf('Duplicate State Found')
            else
%                 sprintf('Inserting New States')
                State=[State,new_state];
                sequence_number = StateTable.insert(new_state);
            end
%             sprintf('Src %d',src)
            labSend({id_parent(2),sequence_number,enabled_transition}, src, 2);
        end
    end
end

for i=1:8
St{i} = State{i};
Ad{i} = Adj{i};
ST{i} = StateTable{i};
end
State = St;
Adj = Ad;
StateTable = ST;
end