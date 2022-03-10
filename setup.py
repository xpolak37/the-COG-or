from distutils.core import setup
from setuptools import setup, find_packages
setup(
  name = 'COGor',
  packages = ['COGor'],
  version = '0.1',
  license='MIT',
  description = 'Package for improving the functional annotation of bacterial genomes, classification '
                'of protein-coding sequences into clusters of orthologous groups, and visualization '
                'of the final annotated genome',
  author = 'Petra Polakovicova',
  author_email = 'xpolak37@vut.cz',
  url = 'https://github.com/xpolak37/the-COG-or',
  download_url = 'https://github.com/xpolak37/the-COG-or/archive/refs/tags/0.1.tar.gz',
  keywords = ['Bacterial genome', 'Functional annotation', 'bioinformatics', 'COG'],
  install_requires=[
          'regex', 'pandas', 'Bio', 'seaborn', 'pillow'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
  ],
  include_package_data=True
)