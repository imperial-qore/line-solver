function java_model = LINE2JLINE(line_model)
switch class(line_model)
    case {'Network', 'MNetwork'}
        java_model = JLINE.from_line_network(line_model);
    case 'LayeredNetwork'
        java_model = JLINE.from_line_layered_network(line_model);
end
end
