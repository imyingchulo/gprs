# setup.py is a script before people lunch your package (in here is "gprs")
# need to specify package name and version, for other information please read:
# https://packaging.python.org/tutorials/packaging-projects/

from setuptools import setup, find_packages

setup(
    name='gprs',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'click>=8.0.0',
        'pandas>=1.1.5'
    ],
    scripts=['bin/gprs']
)