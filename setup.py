import os
import uuid
import glob
from setuptools import setup, find_packages

try: # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError: # for pip <= 9.0.3
    from pip.req import parse_requirements

if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

install_reqs = parse_requirements('requirements.txt', session=uuid.uuid1())
try:
    reqs = [str(ir.requirement) for ir in install_reqs]
except AttributeError:
    reqs = [str(ir.req) for ir in install_reqs]

setup(
    name = "oligotyping",
    version = open('VERSION').read().strip(),
    description = "The oligotyping and minimum entropy decomposition (MED) pipeline for the analysis of marker gene amplicons",
    author = "A. Murat Eren",
    author_email = "meren@mbl.edu",
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
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
    ],

    scripts = [script for script in glob.glob('bin/*') if not script.endswith('-OBSOLETE')],

    include_package_data = True,
    package_data={'': ['Oligotyping/utils/html/scripts/*', 'Oligotyping/utils/html/static/*', 'Oligotyping/utils/html/templates/*']},

    install_requires=reqs,
)
