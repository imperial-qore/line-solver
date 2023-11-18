function varargout = npfqn_traffic_split_cs(MMAP, P, config)
% Given a MMAP, produces a new array after split and class switching
% P(r,(j-1)*R+s): prob that a class-r departure flows to destination j in
% class s out of R possible classes
%empty = cellfun(@isempty, MMAP);
%MMAP(empty)=[];
%P(:,empty)=[];
%P(empty,:)=[];

n = length(MMAP);
[R,J] = size(P);
M = round(J/R);

SMMAP = cell(1,M);
for j=1:M
    SMMAP{j} = cell(1,2+R);
    SMMAP{j}{1} = MMAP{1} + MMAP{2};
    SMMAP{j}{2} = 0*MMAP{2};
    for s=1:R
        SMMAP{j}{2+s} = SMMAP{j}{1} * 0;
        for r=1:R
            SMMAP{j}{2+s} = SMMAP{j}{2+s} + MMAP{2+r}*P(r,(j-1)*R+s);
            SMMAP{j}{2} = SMMAP{j}{2} + MMAP{2+r}*P(r,(j-1)*R+s);
            SMMAP{j}{1} = SMMAP{j}{1} - MMAP{2+r}*P(r,(j-1)*R+s);
        end
    end
    SMMAP{j} = mmap_normalize(SMMAP{j});
end
varargout  = SMMAP;
end