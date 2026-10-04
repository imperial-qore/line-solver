% Test gallery_coxm1 with MVA
model = gallery_coxm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
