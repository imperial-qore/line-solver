function RET = mmap_modulate(P, HT, MMAP)
% MMAP_MODULATE. Modulates in continuous time a set of MMAPs according to
% a given holding time distribution for each MMAP.
% SMMAP = MMAP_MODULATE(HT, MMAP)
% HT{j}: phase-type distribution for holding time in state group j
% MMAPP{j}: departure process while in state group j

J = length(HT);
if length(MMAP)~=J
    error('The HT and MMAP cell arrays must have the same length.');
end

for j=1:length(MMAP)
    if length(MMAP{j})==2 % MAP
        MMAP{j}{3} = MMAP{j}{2}; % autoconvert into MMAP
    end
end

nc = cellfun(@length, MMAP); % nc(j): classes in MMAP{j}
if min(nc) ~= max(nc)
    error('The MMAPs must have the same number of types.');
else
    K = nc(1)-2;
end

nh = cellfun(@(ph) length(ph{1}), HT); % nh(j): states in PH{j}
nm = cellfun(@(mmap) length(mmap{1}), MMAP); % nm(j): states in MMAP{j}

Q = cell(1,J);
RET = cell(1,2+K); % a marked MAP
for k=1:(2+K)
    RET{k} = [];
end
for j=1:J
    Q{j} = cell(1,3+K);
    Q{j}{1} = krons(HT{j}{1}, MMAP{j}{1});
    %Q{j}{2} = krons(HT{j}{2}, MMAP{j}{2});
    Q{j}{3} = kron(HT{j}{2}, eye(nm(j)));
    for k=1:K
        Q{j}{3+k} = kron( eye(nh(j)), MMAP{j}{2+k} );
    end
    nq = length(Q{j}{1});
    MOD_j = cell(1,2+K);
    for k=1:2+K
        MOD_j{k} = [];
    end
    for i=1:J
        if i==j
            MOD_j{1} = [MOD_j{1}, Q{j}{1}];
            for k=1:K
                MOD_j{2+k} = [MOD_j{2+k}, Q{j}{3+k}];
            end            
        else
            entry_i = kron(map_pie(HT{i}),map_pie(MMAP{i}));
            MOD_j{1} = [MOD_j{1}, P(j,i)*Q{j}{3}*ones(nq,1)*entry_i];
            for k=1:K
                MOD_j{2+k} = [MOD_j{2+k}, 0*P(j,i)*Q{j}{3}*ones(nq,1)*entry_i];
            end
        end
    end
    
    RET{1} = [RET{1}; MOD_j{1}];    
    for k=1:K
        RET{2+k} = [RET{2+k}; MOD_j{2+k}];
    end    
end
RET = mmap_normalize(RET);
end