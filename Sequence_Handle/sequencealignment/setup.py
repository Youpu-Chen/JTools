from distutils.core import setup

setup(
  name = 'sequencealignment',
  packages = ['sequencealignment'],
  version = '0.0.4',
  license='MIT',
  description = 'Simple application of Needleman-Wunsch and Smith-Waterman algorithms with affine gap penalty',
  author = 'Youpu Chen',
  author_email = 'otis.hongpu@gmail.com',
  url = 'https://github.com/Youpu-Chen/Myscripts/tree/main/Sequence_Handle/sequencealignment',
  download_url = '',
  keywords = ['bioinformatics', 'sequencealignment', 'Needleman-Wunsch', 'Smith-Waterman', 'affine gap penalty', 'gap extension'],
  install_requires=[
          'numpy',
          'matplotlib',
      ],
  classifiers=[  # Optional
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 3 - Alpha',

    # Indicate who your project is intended for
    'Intended Audience :: Information Technology',
    'Topic :: Scientific/Engineering :: Bio-Informatics',

    # Pick your license as you wish
    'License :: OSI Approved :: MIT License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8',
  ],
)