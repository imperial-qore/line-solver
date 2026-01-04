function classNames = getClassNames(self)
% CLASSNAMES = GETCLASSNAMES()
if self.hasStruct && isfield(self.sn,'classnames')
    classNames = self.sn.classnames;
else
    classNames = string([]); % string array
    for r=1:getNumberOfClasses(self)
        classNames(r,1)=self.classes{r}.name;
    end
end
end