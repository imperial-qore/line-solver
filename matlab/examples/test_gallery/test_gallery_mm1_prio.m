% Test gallery_mm1_prio with MVA
model = gallery_mm1_prio();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
