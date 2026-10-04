function SUP = mmap_super(MMAPa,MMAPb,opt)

% SUP = mmap_super(MMAPa,MMAPb,opt)
if nargin ==2 || (nargin ==3 && strcmpi(opt,'default'))
    % each class in MMAPa and MMAPb is a distinct class in SUP
    K1 = length(MMAPa)-2;
    K2 = length(MMAPb)-2;
    n1 = length(MMAPa{1});
    n2 = length(MMAPb{1});
    SUP = {};
    SUP{1} = krons(MMAPa{1},MMAPb{1});
    SUP{2} = krons(MMAPa{2},MMAPb{2});

    for i=1:K1
        SUP{end+1} = krons(MMAPa{2+i},zeros(n2));
    end

    for j=1:K2
        SUP{end+1} = krons(zeros(n1),MMAPb{2+j});
    end
elseif nargin ==3 && strcmpi(opt,'match')
    % class c in both MMAPa and MMAPb is mapped both into class c of SUP
    K1 = length(MMAPa);
    K2 = length(MMAPb);
    if K1 ~= K2
        error('class matching failed: MMAPs have different number of classes');
    end
    SUP = {};
    for i=1:K1
        SUP{i} = krons(MMAPa{i},MMAPb{i});
    end
elseif nargin == 1
    empty = cellfun(@isempty, MMAPa);
    MMAPa(empty)=[];
    SUP = MMAPa{1};
    for i=2:length(MMAPa)
        SUP = mmap_super(SUP,MMAPa{i});
    end
else
    error('unrecognized option');
end

SUP = mmap_normalize(SUP);
end
