#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages, Extension

laaminoc = Extension(
    name='liblamino',
    extra_compile_args=['-std=c99', '-DWIN32'],
    sources=['lamino.c'])

setup(
    name='lamino',
    packages=find_packages(),
    version=1.0,
    include_package_data=True,
    zip_safe=False,
    author='Doga Gursoy',
    author_email='dgursoy@aps.anl.gov',
    ext_modules=[laaminoc],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: C']
)
