rm dist -rf
python setup.py sdist
python -m twine upload --skip-existing dist/*
