
from setuptools import setup, find_packages

setup(name='cimage_annotation',
      version='2.1.2',
      description='Add cysteine functional annotation to cimage output.',
      author='Dan Bak, Aaron Maurais',
      url='https://github.com/ajmaurais/cimage_annotation',
      classifiers=['Development Status :: 4 - Beta',
        'Intended Audience :: SCIENCE/RESEARCH',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],
      package_dir={'':'src'},
      packages=find_packages(where='src'),
      python_requires='>=3.6.*',
      install_requires=['biopython==1.77.dev0', 'tqdm'],
      dependency_links=['https://github.com/ajmaurais/biopython#egg=biopython-1.77.dev0'],
      entry_points={'console_scripts': ['cimage_annotation=cimage_annotation:main', 'qsub_cimage_annotation=cimage_annotation:qsubmit_main']},
)


