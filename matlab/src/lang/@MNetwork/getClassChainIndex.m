function c = getClassChainIndex(self, jobClass)
% C = GETCLASSCHAININDEX(JOBCLASS)

if ischar(jobClass)
    className = jobClass;
else
    className = jobClass.getName;
end

chains = self.getChains;
for c = 1:length(chains)
    if any(cell2mat(strfind(chains{c}.classnames,className)))
        return
    end
end
c = -1;
end