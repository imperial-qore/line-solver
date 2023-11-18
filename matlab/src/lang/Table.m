function t = Table(varargin)
% T = TABLE(VARARGIN)

% This class was initially introduced for octave compatibility in earlier
% versions
t = table(varargin{:});
for i=1:length(varargin)
    t.Properties.VariableNames{i} = inputname(i);
end
end
