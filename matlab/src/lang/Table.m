function t = Table(varargin)
% T = TABLE(VARARGIN)

t = table(varargin{:});
for i=1:length(varargin)
    t.Properties.VariableNames{i} = inputname(i);
end
end
