from setuptools import setup

setup(name='fitstemplate',
      version='0.1',
      description='Simple package for making FITS template files.',
      packages=['FITSTemplate'],
      install_requires=[
          'numpy',
          'scipy',
          'scikit-image',
          'astropy',
          'future'
      ],
      zip_safe=False)

# setup.py ends here
