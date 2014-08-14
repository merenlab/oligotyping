import os
from setuptools import setup, find_packages
from pip.req import parse_requirements

install_reqs = parse_requirements('requirements.txt')
reqs = [str(ir.req) for ir in install_reqs]

if os.environ.get('USER','') == 'vagrant':
    del os.link

setup(name = "Oligotyping",
      version = "1.0",
      description = "Oligotyping and minimum entropy decomposition (MED) pipeline for the analysis of marker gene amplicons",
      author = u"A. Murat Eren",
      author_email = "meren@mbl.edu",
      url = "http://oligotyping.org",
      packages = find_packages(),
      install_requires=reqs,
      scripts=['decompose', 'entropy-analysis', 'oligotype'],
)

