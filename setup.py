"""
Setup script for helios_postprocess package

Installation:
    pip install -e .  (editable/development mode)
    pip install .     (regular installation)
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "README.md"
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "Helios ICF simulation postprocessing package"

setup(
    name="helios_postprocess",
    version="3.0.0",
    author="Prof T",
    description="ICF post-processing pipeline for Helios EXODUS simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["archive", "archive.*", "examples"]),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "matplotlib>=3.3",
        "netCDF4>=1.5",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "ipython>=7.0",
        ],
        "sklearn": [
            "scikit-learn>=1.0",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords="ICF fusion Helios simulation EXODUS postprocessing",
)
