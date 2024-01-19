#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from setuptools import find_packages, setup

with open("README.md", "r", encoding='utf-8') as fr:
    long_description = fr.read()

required = ['biopython']

setup(
    name="phyloforge",
    version="1.0.0",
    author="Ya Wang",
    author_email="1552082076@qq.com",
    description="PhyloForge: Unifying micro and macro evolution with comprehensive genomic signals",
    license="BSD License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wangyayaya/TreeTool/",
    packages=find_packages(),
    package_data={'': ['*.config', '*.cfg', '*.ini']},
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'phyloforge = treetool.run:main',
        ]
    },
    zip_safe=True,
    install_requires=required
)
