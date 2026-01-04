classdef RoutingMatrix < Copyable
    % Class for routing matrices
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Hidden)
        rt;
    end

    methods
        %Constructor
        function self = RoutingMatrix(par)
            self.rt = par;
        end
        function varargout = subsref(self,s)
            switch s(1).type
                case '.'
                    if length(s) == 1
                        % Implement obj.PropertyName
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    elseif length(s) == 2 && strcmp(s(2).type,'()')
                        % Implement obj.PropertyName(indices)
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    else
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    end
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj(ind).PropertyName
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
                        % Implement obj(indices).PropertyName(indices)
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    else
                        % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    end
                case '{}'
                    if length(s) == 1
                        % Implement obj{indices}
                        [varargout{1:nargout}] = builtin('subsref',self.rt,s);
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj{indices}.PropertyName
                        [varargout{1:nargout}] = builtin('subsref',self.rt,s);
                    else
                        % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',self.rt,s);
                    end
                otherwise
                    error('Not a valid indexing expression')
            end
        end

        function self = subsasgn(self,s,varargin)

            % Allow subscripted assignment to uninitialized variable
            if isequal(self,[])
                % obj = ClassName.empty;
            end

            switch s(1).type
                case '.'
                    if length(s) == 1
                        % Implement obj.PropertyName = varargin{:};
                        self = builtin('subsasgn',self,s,varargin{:});
                    elseif length(s) == 2 && strcmp(s(2).type,'()')
                        % Implement obj.PropertyName(indices) = varargin{:};
                        self = builtin('subsasgn',self,s,varargin{:});
                    else
                        % Call built-in for any other case
                        self = builtin('subsasgn',self,s,varargin{:});
                    end
                case '()'
                    if length(s) == 1
                        % Implement obj(indices) = varargin{:};
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj(indices).PropertyName = varargin{:};
                        self = builtin('subsasgn',self,s,varargin{:});
                    elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
                        % Implement obj(indices).PropertyName(indices) = varargin{:};
                        self = builtin('subsasgn',self,s,varargin{:});
                    else
                        % Use built-in for any other expression
                        self = builtin('subsasgn',self,s,varargin{:});
                    end
                case '{}'
                    if length(s) == 1
                        % Implement obj{indices} = varargin{:}
                        if isa(varargin{1},'jline.lang.RoutingMatrix')
                            jobclass = varargin{1}.getClass(s.subs{1}-1);
                            mat = JLINE.from_jline_matrix(varargin{1}.get(jobclass,jobclass));
                            self.rt = builtin('subsasgn',self.rt,s,mat);
                        else
                            self.rt = builtin('subsasgn',self.rt,s,varargin{1});
                        end
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj{indices}.PropertyName = varargin{:}
                        self = builtin('subsasgn',self.rt,s,varargin{:});
                    else
                        % Use built-in for any other expression
                        self = builtin('subsasgn',self.rt,s,varargin{:});
                    end
                otherwise
                    error('Not a valid indexing expression')
            end
        end

        function self = set(self, varargin)
            jobclass1=varargin{1};
            jobclass2=varargin{2};
            if length(varargin)==5
                node1=varargin{3};
                node2=varargin{4};
                val=varargin{5};
                self.rt{jobclass1,jobclass2}(node1, node2)=val;
            elseif length(varargin)==3
                mat=varargin{3};
                self.rt{jobclass1,jobclass2}=mat;
            elseif length(varargin)==2
                mat=varargin{2};
                self.rt{jobclass1,jobclass1}=mat;
            else
                line_error(mfilename,'Invalid number of arguments.');
            end
        end

        function P = getCell(self)
            P = self.rt;
        end

        function print(self)
            for r = 1:size(self.rt,1)
                for s = 1:size(self.rt,2)
                    fprintf(1,'Class %d -> Class %d:\n',r,s);
                    disp(self.rt{r,s});
                end
            end
        end
    end

    methods (Static)
        function [rtorigcell,rtorig] = rtnodes2rtorig(sn)
            K = sn.nclasses;
            rtnodes = sn.rtnodes;

            csshift = sn.nnodes;
            for ind = 1:sn.nnodes
                if startsWith(sn.nodenames{ind}, 'CS_')
                    csshift = ind-1;
                    break
                end
            end

            colToKeep=[];
            for ind = 1:csshift
                for k = 1:K
                    colToKeep(end+1) = (ind-1)*K+k;
                end
            end

            rtorig = dtmc_stochcomp(rtnodes, colToKeep);
            rtorigcell = cellzeros(K,K,csshift,csshift);

            % for cache rt, replace NaNs for unknown probabilities with 0
            rtorig(isnan(rtorig)) = 0;

            for ind = 1:csshift
                if sn.nodetype(ind) ~= NodeType.Sink
                    for jnd = 1:csshift
                        for r = 1:K
                            for s = 1:K
                                rtorigcell{r,s}(ind,jnd) = rtorig((ind-1)*K+r,(jnd-1)*K+s);
                            end
                        end
                    end
                end
            end

        end
    end
end