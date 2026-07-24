function n = get_max_factor(varargin)
if nargin == 0
    n = 177;
    return;
end

type_name = varargin{1};
if isa(type_name, 'string') || isa(type_name, 'char')
    switch char(type_name)
        case 'single'
            n = 38;
        case 'double'
            n = 177;
        otherwise
            n = 255;
    end
else
    if isa(type_name, 'single')
        n = 38;
    elseif isa(type_name, 'double')
        n = 177;
    else
        n = 255;
    end
end
end
