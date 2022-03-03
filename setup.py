#!/usr/bin/env python

#version: 1.0
#author: Evgeniy Mozheiko
#Contact: eamozheiko@gmail.com

import os,sys,io
from distutils.core import setup, Extension
from setuptools import find_packages
from os import path


this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


requires = []

def main():
    #compilemis()
    setup(name="ExoC-Package",
          version="1.0",
          description="ExoC - Translocation detection from Exome capture Hi-C",
          long_description=long_description,
          long_description_content_type='text/markdown',
          author='Evgeniy Mozheiko',
          author_email='eamozheiko@gmail.com',
          package_dir={'ExoC' : 'ExoC'},
          install_requires = requires,
          setup_requires = requires,
          #packages=['HiNT'],
          packages=find_packages(),
          package_data={'ExoC':['scripts/*']},
          include_package_data=True,
          scripts=['bin/exoc'],
	  url="https://github.com/eamozheiko/ExoC",
          classifiers=[
            'Programming Language :: Python :: 3',
            'Environment :: Console',
            'Environment :: Web Environment',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: GPL-3.0 License',
            'Operating System :: Linux',
            'Topic :: Software Development',
            ],
          )


if __name__ == '__main__':
    main()
