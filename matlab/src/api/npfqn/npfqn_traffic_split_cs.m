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
for jst=1:M
    SMMAP{jst} = cell(1,2+R);
    SMMAP{jst}{1} = MMAP{1} + MMAP{2};
    SMMAP{jst}{2} = 0*MMAP{2};
    for s=1:R
        SMMAP{jst}{2+s} = SMMAP{jst}{1} * 0;
        for r=1:R
            SMMAP{jst}{2+s} = SMMAP{jst}{2+s} + MMAP{2+r}*P(r,(jst-1)*R+s);
            SMMAP{jst}{2} = SMMAP{jst}{2} + MMAP{2+r}*P(r,(jst-1)*R+s);
            SMMAP{jst}{1} = SMMAP{jst}{1} - MMAP{2+r}*P(r,(jst-1)*R+s);
        end
    end
    SMMAP{jst} = mmap_normalize(SMMAP{jst});
end
varargout  = SMMAP;
end