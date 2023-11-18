function options=Q_ParseOptPara(options,Names,Types,Values,OptParas)

for i=1:2:size(OptParas,2)
    prop=OptParas{i};
    if (i+1 <= size(OptParas,2)) 
        value=OptParas{i+1};
        arg=strmatch(prop,Names,'exact');
        if (isempty(arg)) % checks whether Invalid Property
             warning('MATLAB:Q_ParseOptPara:InvalidPropName',...
                'Property name ''%s'' not recognized and ignored',prop);
        else  
            if (eval(['is' Types(arg,:) '(value)'])) % checks whether property value is of correct type
                if (strmatch('char',Types(arg,:),'exact')) % property values are strings
                    %if (isempty(strmatch(value,eval([prop 'Value']),'exact')))
                    if (isempty(strmatch(value,Values{arg},'exact')))
                        warning('MATLAB:Q_ParseOptPara:InvalidPropValue',...
                        'Property value ''%s'' of ''%s'' not allowed and ignored',value,prop);
                    else    
                        options.(prop)=value;
                    end
                elseif (strmatch('numeric',Types(arg,:),'exact'))   % property values are numeric
                    options.(prop)=value;
                elseif (strmatch('cell',Types(arg,:),'exact'))   % property values are cell
                    options.(prop)=value;
                end
            else % incorrect property value type
                if (ischar(value))
                    warning('MATLAB:Q_ParseOptPara:InvalidPropType',...
                    'Property value ''%s'' of ''%s'' has an incorrect type and is ignored',value,prop);
                end
                if (iscell(value))
                    warning('MATLAB:Q_ParseOptPara:InvalidPropType',...
                    'Property value ''%s'' of ''%s'' has an incorrect type and is ignored',value,prop);
                end
                if (isnumeric(value))
                    warning('MATLAB:Q_ParseOptPara:InvalidPropType',...
                    'Property value %d of ''%s'' has an incorrect type and is ignored',value,prop);
                end
            end    
        end
    else % odd number of optional parameters
        warning('MATLAB:Q_ParseOptPara:OddNumbOptParas',...
            'An odd number of optional parameters detected, last parameter ''%s'' ignored',prop);
    end    
end
