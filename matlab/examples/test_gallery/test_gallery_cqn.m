% Test gallery_cqn with MVA
model = gallery_cqn();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
