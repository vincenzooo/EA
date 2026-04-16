from setuptools import setup

# see https://docs.python.org/2/distutils/examples.html#pure-python-distribution-by-package
setup(
  name='EA',
  version='0.0.1',
  description="Calculate Reflectivity and Effective area of a nested shell X-ray telescope.",
  author='Vincenzo Cotroneo',
  author_email='vincenzo.cotroneo@inaf.it',
  install_requires=['wheel','numpy','matplotlib','scipy','IPython'],
  package_dir={'': '.'},
  setup_requires=['numpy'],
  include_package_data=True
)

