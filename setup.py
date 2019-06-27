#!/usr/bin/env python3

import setuptools

with open("README.md",'r') as f:
    long_description = f.read()

setuptools.setup(
        name="martini-garnish",
        version="0.2-alpha",
        author=["Matthijs Tadema", "Lorenzo Gaifas"],
        author_email="M.J.Tadema@protonmail.com",
        description="Render coarse grained molecular structures in PyMOL",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/mjtadema/pycg_bonds",
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            ],
        install_requires=[
            'networkx',
            'numpy'
            ],
        python_requires='>=3.5'
        )
