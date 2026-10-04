% Test gallery_parm1 with MVA
model = gallery_parm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
