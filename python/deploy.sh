rm dist -rf
python setup.py sdist
python -m twine upload dist/*
