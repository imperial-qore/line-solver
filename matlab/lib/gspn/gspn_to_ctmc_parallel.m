function [State, Adj, StateTable] = gspn_to_ctmc_parallel(F, B, INH, IM, P, M0)


[num_m, num_t] = size(F);

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
    C = cellfun(@minus,B,F,'Un',0);
    i=0;
    
    while e < 10000
        e = e + 1;
        if i < size(State,2)
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
                assigned_node = h0(new_state, N);
                if assigned_node == labindex
                    [exists,sequence_number] = StateTable.search(new_state);
                    if exists
                        sprintf('Duplicate State Found')
                        Adj.insert(i, [labindex, sequence_number], [enabled_m(j),enabled_t(j)]);
                    else
                        State=[State,new_state];
                        sequence_number = StateTable.insert(new_state);
                        Adj.insert(i, [labindex, sequence_number], [enabled_m(j),enabled_t(j)]);
                    end
                else
                    fprintf('Sent Data %s\n',class(new_state))
                    labSend({[labindex, i],new_state, [enabled_m(j),enabled_t(j)]}, assigned_node, 1); % Tag == 1 Send State
                end
            end
        end
        while labProbe('any', 2) % Tag == 2 IDs
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