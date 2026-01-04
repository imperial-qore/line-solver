function P = getLinkedRoutingMatrix(self)
% P = GETLINKEDROUTINGMATRIX()

sn = self.getStruct(false);
P = sn.rtorig;

% if isempty(sn.rtorig)
%     %line_warning(mfilename,'Unsupported: getLinkedRoutingMatrix() requires that the model topology has been instantiated with the link() method. Attempting auto-recovery.\n');
%     %fname = lineTempName;
%     %QN2MATLAB(self,self.name,fname);
%     %run(fname);
%     %delete(fname);
%     %P = model.getStruct().rtorig;
% end

end