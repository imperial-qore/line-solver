% Test gallery_mhyp1 with MVA
model = gallery_mhyp1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
