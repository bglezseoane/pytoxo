# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python tool to calculate penetrance tables for
# high-order epistasis models
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

import setuptools


with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

with open("requirements_with_gui.txt") as f:
    requirements = f.read().splitlines()


version = "1.0"

setuptools.setup(
    name="pytoxo",
    version=version,
    packages=["pytoxo", "pytoxo_cli", "pytoxo_gui"],
    entry_points={
        "console_scripts": [
            "pytoxo_cli = pytoxo_cli.__main__:main",
            "pytoxo_gui = pytoxo_gui.__main__:main",
        ],
    },
    python_requires=">=3.8",
    install_requires=requirements,
    data_files=[
        ("", ["requirements_with_gui.txt"]),
        ("", ["README.md"]),
        ("", ["LICENSE"]),
    ],
    url="https://github.com/bglezseoane/pytoxo",
    download_url=f"https://github.com/bglezseoane/pytoxo/archive/{version}.tar.gz",
    license="LICENSE",
    author="Borja González Seoane",
    author_email="borja.gseoane@udc.es",
    description="A Python tool to calculate penetrance tables for high-order epistasis models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
)
