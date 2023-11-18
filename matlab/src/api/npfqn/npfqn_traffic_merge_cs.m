function SMMAP = npfqn_traffic_merge_cs(MMAP, prob, config)
% Given a cell array of n MMAPs with arrivals of R classes, produces a new
% MMAP after flow merging  and class switching
% prob((i-1)*R+r,s): prob that a class-r arrival from MMAP i switches to class-s
n = length(MMAP);
[~,R] = size(prob);

if ~isfield(config,'merge')
    config.merge = 'default';
end

% First let's convert from types to classes
for i=1:n
    P = zeros(R);
    for r=1:R
        for s=1:R
            P(r,s) = prob((i-1)*R+r,s);
        end
    end
    MMAP{i} = mmap_mark(MMAP{i}, P);
end

if n == 1
    SMMAP = MMAP{1};
    return;
else
    switch config.merge
        case {'default','super'}
            SMMAP = MMAP{1};
            for j=2:n
                SMMAP = mmap_super(SMMAP, MMAP{j},'match');
            end
    end
end
end