import re
import setuptools
import subprocess
import sys
import os

from biosimulators_smoldyn import __version__


# set dirs 
name = 'biosimulators_smoldyn'
dirname = os.path.dirname(__file__)


# set descriptions
with open("README.md", "r") as readme:
    description = readme.read()
    # Patch the relative links to absolute URLs that will work on PyPI.
    description2 = re.sub(
        r']\(([\w/.-]+\.png)\)',
        r'](https://github.com/biosimulators/Biosimulators_smoldyn/raw/main/\1)',
        description)
    long_description = re.sub(
        r']\(([\w/.-]+)\)',
        r'](https://github.com/biosimulators/Biosimulators_smoldyn/blob/main/\1)',
        description2)


# install package
setuptools.setup(
    name=name,
    version=__version__,
    description=("BioSimulators-compliant command-line interface to the Smoldyn simulation program <https://github.com/ssandrews/Smoldyn>."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/biosimulators/Biosimulators_Smoldynå",
    download_url="https://github.com/biosimulators/Biosimulators_Smoldyn",
    author='Center for Reproducible Biomedical Modeling',
    author_email="info@biosimulators.org",
    license="MIT",
    keywords=['BioSimulators', 'systems biology', 'computational biology', 'mathematical model',
              'kinetic model', 'simulation', 'stochastic', 'spatial', 'SED-ML', 'COMBINE', 'OMEX'],
    packages=setuptools.find_packages(exclude=['tests', 'tests.*']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    install_requires=[
        "smoldyn>=2.73",
        "biosimulators-utils[logging]>=0.1.188",
        "simulariumio>=1.11.0",
        "pkg-utils"
    ],
    entry_points={
        'console_scripts': [
            'biosimulators-amici = biosimulators_smoldyn.__main__:main',
        ],
    },
)
