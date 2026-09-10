function lbl=label(varargin)
if length(varargin)==1
    cellchar = varargin{1};
    lbl=categorical(cellchar);
else
    n=varargin{1};
    m=varargin{2};
    lbl=categorical(n,m);
end
end