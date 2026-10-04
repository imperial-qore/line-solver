% Test gallery_repairmen with MVA
model = gallery_repairmen();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
