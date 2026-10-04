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
            % A Join that names no Fork is REFUSED rather than left as an
            % empty row. The pairing is a declaration, not a derivation from
            % the routing (examples/basic/forkJoin/fj_basic_nesting has two
            % forks and two joins whose pairing the routing alone does not
            % fix), so an unpaired Join never fires and every consumer of sn.fj
            % would go on to analyze a model the user did not write. This used
            % to die one line below on "[] has no property index"; the JAR
            % warned and returned the empty row and native python was silent.
            if isempty(self.nodes{ind}.joinOf)
                line_error(mfilename, sprintf(['Join ''%s'' closes no Fork: pass the Fork ' ...
                    'to the Join constructor, Join(model, name, fork).'], self.nodes{ind}.name));
            end
            fjPairs(self.nodes{ind}.joinOf.index,self.nodes{ind}.index) = true;
    end
end
ft = fjPairs;
end
