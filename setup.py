#!/usr/bin/env python3
import setuptools

with open("README.md", 'r') as f:
    long_description = f.read()

setuptools.setup(
        name="garnish",
        version="0.4-alpha",
        author=["Matthijs Tadema", "Lorenzo Gaifas"],
        author_email=["M.J.Tadema@protonmail.com", "brisvag@gmail.com"],
        description="Render coarse grained molecular structures in PyMOL",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/mjtadema/garnish",
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            ],
        python_requires='>=3.5',
        package_data={
            "garnish": ["data/*"]
            }
        )
