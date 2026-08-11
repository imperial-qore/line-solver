% Test gallery_mapm1 with MVA
model = gallery_mapm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
