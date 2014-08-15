import os
import glob
from setuptools import setup, find_packages
from pip.req import parse_requirements


if os.environ.get('USER','') == 'vagrant':
    del os.link

os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))

with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()

install_reqs = parse_requirements('requirements.txt')
reqs = [str(ir.req) for ir in install_reqs]

setup(
    name = "Oligotyping",
    version = "1.0",
    description = "Oligotyping and minimum entropy decomposition (MED) pipeline for the analysis of marker gene amplicons",
    author = u"A. Murat Eren",
    author_email = "meren@mbl.edu",
    license = "GPLv3+",
    url = "http://oligotyping.org",
    packages = find_packages(),

    longer_description=README,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: R',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
    ],

    scripts = [script for script in glob.glob('bin/*') if not script.endswith('-OBSOLETE')],

    install_requires=reqs,
)
