rm dist/* -rf
rm line_solver.egg_info -rf
python setup.py sdist
twine upload dist/*
