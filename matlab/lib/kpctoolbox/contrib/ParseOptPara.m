function options=ParseOptPara(options,Names,Types,Values,OptParas)
% Script distributed with "D. A. Bini, B. Meini, S. Steffe, B. Van Houdt: 
% Structured Markov chains solver: software tools"

for i=1:2:size(OptParas,2)
    prop=OptParas{i};
    if (i+1 <= size(OptParas,2)) 
        value=OptParas{i+1};
        arg=find(strcmp(prop,cellstr(Names)));
        if (isempty(arg)) % checks whether Invalid Property
             warning('MATLAB:ParseOptPara:InvalidPropName',...
                'Property name ''%s'' not recognized and ignored',prop);
        else  
            if (eval(['is' Types(arg,:) '(value)'])) % checks whether property value is of correct type
                if (strcmp('char',deblank(Types(arg,:)))) % property values are strings
                    %if (isempty(strmatch(value,eval([prop 'Value']),'exact')))
                    if (~any(strcmp(value,cellstr(Values{arg}))))
                        warning('MATLAB:ParseOptPara:InvalidPropValue',...
                        'Property value ''%s'' of ''%s'' not allowed and ignored',value,prop);
                    else    
                        options.(prop)=value;
                    end
                elseif (strcmp('numeric',deblank(Types(arg,:))))   % property values are numeric
                    options.(prop)=value;
                end
            else % incorrect property value type
                if (ischar(value))
                    warning('MATLAB:ParseOptPara:InvalidPropType',...
                    'Property value ''%s'' of ''%s'' has an incorrect type and is ignored',value,prop);
            
                end
                if (isnumeric(value))
                    warning('MATLAB:ParseOptPara:InvalidPropType',...
                    'Property value %d of ''%s'' has an incorrect type and is ignored',value,prop);
            
                end
            end    
        end
    else % odd number of optional parameters
        warning('MATLAB:ParseOptPara:OddNumbOptParas',...
            'An odd number of optional parameters detected, last parameter ''%s'' ignored',prop);
    end    
end
