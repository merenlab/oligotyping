import os
import glob
from setuptools import setup, find_packages

if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

# read the dependencies straight from requirements.txt. we deliberately avoid importing
# anything from `pip` here: pip's internal API is unsupported and moves around between
# releases, and it is not even importable inside the isolated build environment that
# modern pip uses when building the package.
with open('requirements.txt') as f:
    reqs = [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name = "oligotyping",
    version = open('VERSION').read().strip(),
    description = "The oligotyping and minimum entropy decomposition (MED) pipeline for the analysis of marker gene amplicons",
    author = "A. Murat Eren",
    author_email = "a.murat.eren@gmail.com",
    license = "GPLv3+",
    url = "http://oligotyping.org",
    packages = find_packages(),

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],

    python_requires = '>=3.7',

    scripts = [script for script in glob.glob('bin/*') if not script.endswith('-OBSOLETE')],

    include_package_data = True,
    package_data={'': ['Oligotyping/utils/html/scripts/*', 'Oligotyping/utils/html/static/*', 'Oligotyping/utils/html/templates/*']},

    install_requires=reqs,
)
