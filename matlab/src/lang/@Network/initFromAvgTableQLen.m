function initFromAvgTableQLen(self, AvgTable)
QN = reshape(AvgTable.QLen,self.getNumberOfClasses,self.getNumberOfStations)';
self.initFromAvgQLen(QN);
end