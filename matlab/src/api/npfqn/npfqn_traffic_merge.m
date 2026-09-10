function SMMAP = npfqn_traffic_merge(MMAP, config)
%Given a cell array of n MMAPs with arrivals of R classes, produces a new
%MMAP after flow merging
empty = cellfun(@isempty, MMAP);
MMAP(empty)=[];
n = length(MMAP);

if nargin<2%~exist('config','var')
    config = struct('merge','default');
end

if n == 1
    SMMAP = MMAP{1};
    return;
else
    switch config.merge
        case {'default','super'}
            SMMAP = MMAP{1};
            for j=2:n
                if ~isempty(MMAP{j})
                    SMMAP = mmap_super(SMMAP, MMAP{j},'match');
                end
            end
        case 'mixture'
            SMMAP = MMAP{1};
            for j=2:n
                if ~isempty(MMAP{j})
                    SMMAP = mmap_super(SMMAP, MMAP{j},'match');
                    if mmap_isfeasible(SMMAP)
                        SMMAP = mmap_mixture_fit_mmap(SMMAP);
                    end
                end
            end
        case 'interpos'
            FLOW = {};
            for j=1:n
                if ~isempty(MMAP{j})
                    FLOW{end+1} = m3pp2m_fitc_theoretical(MMAP{j}, 'exact_delta', 1, 1e6);
                end
            end
            SMMAP = m3pp2m_interleave(FLOW);
        otherwise
            line_error(mfilename,'Unsupported configuration for merge.');
    end
    
    switch config.compress
        case 'none'
            % do nothing
        case 'default'
            SMMAP = mmap_compress(SMMAP);
    end
end
SMMAP = mmap_normalize(SMMAP);
end
