function java_model = LINE2JLINE(line_model)
% JAVA_MODEL = LINE2JLINE(LINE_MODEL) the jline mirror of a Network, JNetwork or LayeredNetwork
%
% JLINE.line_to_jline forwards here, so this is the single implementation.
if isa(line_model, 'MNetwork')
    java_model = JLINE.from_line_network(line_model);
elseif isa(line_model, 'JNetwork')
    java_model = line_model.obj;
elseif isa(line_model, 'LayeredNetwork')
    java_model = JLINE.from_line_layered_network(line_model);
else
    line_error(mfilename, sprintf('LINE2JLINE supports Network, JNetwork and LayeredNetwork, got a %s.', class(line_model)));
end
end
