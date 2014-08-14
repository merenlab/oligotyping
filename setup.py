import os
from setuptools import setup, find_packages
from pip.req import parse_requirements

install_reqs = parse_requirements('requirements.txt')
reqs = [str(ir.req) for ir in install_reqs]

if os.environ.get('USER','') == 'vagrant':
    del os.link

setup(name = "Oligotyping",
      version = "0.1",
      description = "Oligotyping library for Python",
      author = u"A. Murat Eren",
      author_email = "meren@mbl.edu",
      url = "https://github.com/meren/oligotyping",
      packages = find_packages(),
      install_requires=reqs,
      scripts=['decompose', 'entropy-analysis', 'oligotype'],
)

