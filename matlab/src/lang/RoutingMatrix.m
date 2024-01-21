classdef RoutingMatrix < Copyable
    % Class for routing matrices
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties (Hidden)
        P;
    end

    methods
        %Constructor
        function self = RoutingMatrix(Pcell)
            self.P = Pcell;
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
                        [varargout{1:nargout}] = builtin('subsref',self.P,s);
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj{indices}.PropertyName
                        [varargout{1:nargout}] = builtin('subsref',self.P,s);
                    else
                        % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',self.P,s);
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
                        self = builtin('subsasgn',self.P,s,varargin{:});
                    elseif length(s) == 2 && strcmp(s(2).type,'.')
                        % Implement obj{indices}.PropertyName = varargin{:}
                        self = builtin('subsasgn',self.P,s,varargin{:});
                    else
                        % Use built-in for any other expression
                        self = builtin('subsasgn',self.P,s,varargin{:});
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
                self.P{jobclass1,jobclass2}(node1, node2)=val;
            else
                mat=varargin{3};
                self.P{jobclass1,jobclass2}=mat;
            end
        end

        function P = getCell(self)
            P = self.P;
        end
    end
end