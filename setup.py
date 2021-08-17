
from setuptools import setup, find_packages

##def readme():
##    with open('README.rst') as f:
##        return f.read()

setup(name='krige_tools',
      version='0.1',
      description='krige tools',
##      long_description=readme(),
##      url='https://pypi.python.org/pypi/zbridge',
      author='Matthias Preisig',
      author_email='mpreisig@geomod.ch',
      license='Geomod',
      classifiers=['Topic :: Scientific/Engineering'],
      packages=find_packages(exclude=['dist']),
      install_requires=[
      ],
      include_package_data=True,
      zip_safe=False)
