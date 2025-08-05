#!/usr/bin/env python3
"""
Setup script for the dislocate package.
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements from requirements.txt
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        requirements = []
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                requirements.append(line)
        return requirements

setup(
    name="dislocate",
    version="0.1.0",
    author="David Jany",
    author_email="djany31@gmail.com",
    description="Python tools for analyzing dislocation core structures",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/dj-shark31/dislocate",  # Update with actual repository URL
    project_urls={
        "Bug Tracker": "https://github.com/dj-shark31/dislocate/issues",  # Update with actual URL
        "Documentation": "https://github.com/dj-shark31/dislocate#readme",  # Update with actual URL
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Materials Science",
        "License :: OSI Approved :: MIT License",  # Update based on your LICENSE file
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="dislocation, materials science, atomistic simulation, crystal defects, physics",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
        ],
        "lammps": [
            "lammps>=2022.0",
        ],
    },
    entry_points={
        "console_scripts": [
            # Add any command-line tools here if you have them
            # "dislocate=dislocate.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "dislocate": [
            "config.example.yaml",
            "sqs_cells/**/*.out",
            "sqs_cells/**/*.in",
            "potentials/*",
            "examples/**/*",
        ],
    },
    exclude_package_data={
        "": [
            "*.pyc",
            "__pycache__",
            ".DS_Store",
            ".git*",
            "test/*",
        ],
    },
    zip_safe=False,
) 