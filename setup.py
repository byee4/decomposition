import os
from setuptools import setup


# Utility function to read the README.md file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README.md file and 2) it's easier to type in the README.md file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "decompose",
    version = "0.0.1",
    author = "Brian",
    author_email="byee4@jhu.edu",
    description = ("decomposition tools"),
    license = "BSD",
    keywords = "decomposition",
    url = "http://github.com/byee4/decomposition",
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
    packages=find_packages(), # ['decomposition'],
    package_dir={'decomposition': 'decomposition/'},
    include_package_data=True,
    install_requires=[
        'numpy>=1.10',
        'pandas>=0.16',
        'matplotlib>=1.5',
        'sklearn',
        'seaborn>=0.7',
        'bokeh>=0.10.0'
    ],
    entry_points = {
        'console_scripts': [
            'decompose = decomposition.decompose:main'
        ]
    }
)