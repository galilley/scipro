#!/usr/bin/env python

metadata = dict(
    name="scipro",
    version="0.12",
    license="GNU GPL",
    platforms="POSIX",
    author="Denis Kharenko",
    author_email="kharenko@iae.nsk.su",
    description=("Scientific Processing toolkit"),
    keywords=["photonic", "fiber", "fibre", "optical", "equipment"],
    packages=["scipro", "scipro/reader"],
    long_description=open("README.md").read(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics"
    ]
)

if __name__ == "__main__":

    try:
        from setuptools import setup

    except ImportError:
        from distutils.core import setup

    # No matter which imports were used, use same setup call:
    setup(**metadata)
