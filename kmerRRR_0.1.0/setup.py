#!/usr/bin/env python3

from setuptools import setup, find_packages #type: ignore

setup(
    name="kmerRRR",
    version="0.1.0",
    description="k-mer based tool for functional genomics in Repeat Rich Regions",
    author = "Jabale Rahmat",
    author_email="jrahmat@ur.rochester.edu",
    packages=find_packages(), include_package_data= True, package_data= {"kmerRRR": ["R/pbr_plot.r", "R/mapq_plot.r"]},
    install_requires=["numpy", "biopython", "scipy", "pysam", "matplotlib", "pyjellyfish"],
    entry_points={
        'console_scripts': [
            'kmerRRR = kmerRRR.cli:main',
        ],
    },
    python_requires='>=3.6'
)