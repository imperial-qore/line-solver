function Q_CheckUnusedParaMG1(options)
%Q_CheckUnusedPara checks if the optional parameter include unused entries
if ( isempty(strfind(options.Mode,'CR')) )&&( ~isempty(options.OptMG1_CR) )
    warning('MATLAB:Q_CheckUnusedParaMG1:IgnoredParamateter',...
    'Optional parameter ''OptMG1_CR'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'FI')) )&&( ~isempty(options.OptMG1_FI) )
    warning('MATLAB:Q_CheckUnusedParaMG1:IgnoredParamateter',...
    'Optional parameter ''OptMG1_FI'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'IS')) )&&( ~isempty(options.OptMG1_IS) )
    warning('MATLAB:Q_CheckUnusedParaMG1:IgnoredParamateter',...
    'Optional parameter ''OptMG1_IS'' ignored under mode ''%s''',options.Mode);
elseif( isempty(strfind(options.Mode,'RR')) )&&( ~isempty(options.OptMG1_RR) )
    warning('MATLAB:Q_CheckUnusedParaMG1:IgnoredParamateter',...
    'Optional parameter ''OptMG1_RR'' ignored under mode ''%s''',options.Mode);
end