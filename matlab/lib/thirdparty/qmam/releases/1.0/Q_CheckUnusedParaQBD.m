function Q_CheckUnusedParaQBD(options)
%Q_CheckUnusedPara checks if the optional parameter include unused entries
if ( isempty(strfind(options.Mode,'CR')) )&&( ~isempty(options.OptQBD_CR) )
    warning('MATLAB:Q_CheckUnusedParaQBD:IgnoredParamateter',...
    'Optional parameter ''OptQBD_CR'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'FI')) )&&( ~isempty(options.OptQBD_FI) )
    warning('MATLAB:Q_CheckUnusedParaQBD:IgnoredParamateter',...
    'Optional parameter ''OptQBD_FI'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'IS')) )&&( ~isempty(options.OptQBD_IS) )
    warning('MATLAB:Q_CheckUnusedParaQBD:IgnoredParamateter',...
    'Optional parameter ''OptQBD_IS'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'LR')) )&&( ~isempty(options.OptQBD_LR) )
    warning('MATLAB:Q_CheckUnusedParaQBD:IgnoredParamateter',...
    'Optional parameter ''OptQBD_LR'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'NI')) )&&( ~isempty(options.OptQBD_NI) )
    warning('MATLAB:Q_CheckUnusedParaQBD:IgnoredParamateter',...
    'Optional parameter ''OptQBD_NI'' ignored under mode ''%s''',options.Mode);
end
