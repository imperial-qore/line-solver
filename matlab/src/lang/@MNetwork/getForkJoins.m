function ft = getForkJoins(self)
% FT = GETFORKJOINS()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

I = getNumberOfNodes(self);

%K = getNumberOfClasses(self);
% ft = zeros(M*K); % fork table
% for i=1:M % source
%     for r=1:K % source class
%         for j=1:M % dest
%             switch class(self.stations{i})
%                 case 'Fork'
%                     if rt((i-1)*K+r,(j-1)*K+r) > 0
%                         ft((i-1)*K+r,(j-1)*K+r) = self.stations{i}.output.tasksPerLink;
%                     end
%             end
%         end
%     end
% end

fjPairs = false(I,I);
for ind=1:I
    switch class(self.nodes{ind})
        case 'Fork'
            % no-op
        case 'Join'
            fjPairs(self.nodes{ind}.joinOf.index,self.nodes{ind}.index) = true;
    end
end
ft = fjPairs;
end
