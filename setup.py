"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
from xmipp2 import __version__

here = path.abspath(path.dirname(__file__))

# Load requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()


# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-xmipp2',  # Required
    version=__version__,  # Required
    description='xmipp2 ready to use in scipion.',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/I2PC/scipion-em-xmipp2.git',  # Optional
    author='I2PC',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion cryoem imageprocessing scipion-3.0',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={  # Optional
       'pyworkflow.plugin': 'xmipp2 = xmipp2',
    },
    package_data={  # Optional
       'xmipp2': ['xmipp_logo.png', 'protocols.conf'],
    }
)
