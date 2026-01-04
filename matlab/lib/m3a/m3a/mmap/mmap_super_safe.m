function SUP = mmap_super_safe(MMAPs, maxorder, method)
% SUP = mmap_super_safe(MMAPs, maxorder)
% Superposition of marked MAPS

if nargin == 2
    method = 'default';
end
empty = cellfun(@isempty, MMAPs);
MMAPs(empty)=[];
scv_unmarked = cellfun(@map_scv, MMAPs);
[~,Iset]=sort(scv_unmarked); % sort flows with small scv first
SUP = {};
for i=Iset(:)' % low-SCV first
    if isempty(SUP) % is this is the first MMAP
        SUP = MMAPs{i};
        if maxorder == 1
            %  then treat it as a Poisson process if the limit is 1
            SUP = mmap_exponential(mmap_lambda(MMAPs{i}));
        end
    else
        if length(SUP{1}) * length(MMAPs{i}{1}) > maxorder
            %  otherwise treat it as a marked AMAP(2) process if the limit is >1
            if length(SUP{1}) * 2 <= maxorder
                SUP = mmap_super(SUP, mamap2m_fit_gamma_fb_mmap(MMAPs{i}), method);
            else
                SUP = mmap_super(SUP, mmap_exponential(mmap_lambda(MMAPs{i})), method);
            end
        else
            SUP = mmap_super(SUP,MMAPs{i}, method);
        end
    end
  
end
end