from setuptools import setup, find_packages

setup(
    name='line_solver',
    version='2.0.31.10',
    author='Giuliano Casale',
    author_email='g.casale@imperial.ac.uk',
    description='Queueing theory algorithms in Python',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
