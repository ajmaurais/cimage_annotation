
from distutils.core import setup

setup(name='cimage_annotation',
      version='2.0.0',
      description='Add cysteine functional annotation to cimage output.',
      author='Dan Bak, Aaron Maurais',
      url='https://github.com/ajmaurais/cimage_annotation',
      classifiers=['Development Status :: 3 - Alpha',
        'Intended Audience :: SCIENCE/RESEARCH',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],
      package_dir={'':'src'},
      packages=find_packages(where='src'),
      python_requires='>=3.6.*, >=3.7.*, >=3.8.*',
      install_requires=['biopython'],
      dependency_links=['https://github.com/ajmaurais/biopython/tarball/swissprot_bugfix#egg=package-1.0'],
      entry_points={'console_scripts': ['cimage_annotation=cimage_annotation:main']}
)


