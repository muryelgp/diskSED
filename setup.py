from setuptools import setup, find_packages
import os.path

def read(rel_path):
  here = os.path.abspath(os.path.dirname(__file__))
  with open(os.path.join(here, rel_path), 'r', encoding='utf8') as fp:
    return fp.read()


setup(
    name='diskSED',
    version='1.0.0',
    description='Model to fit full SED (X-ray spectra and UV/opt data) of compact accretion disk around black holes. Includes also other tools for fully exploration of the fitting results.',
    long_description=read('README.md'),
    author='Muryel Guolo',
    author_email='mguolop1@jhu.edu',
    url='https://github.com/muryelgp/diskSED',
    download_url = 'https://github.com/muryelgp/diskSED.git',
    packages=['diskSED'],
    install_requires=[
          'numpy',
          'matplotlib',
	'astropy',
      ],
)