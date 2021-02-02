#!/usr/bin/env python

# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python library for calculating penetrance tables of any
# bivariate epistasis model.
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

"""Run all unit tests helper script.

Simple script to discover and run all tests inside the current directory
with the project home directory as working directory."""

import os
import unittest

# Assert project home directory because tests have routes since here
if os.path.basename(os.getcwd()) == "test":
    os.chdir(os.pardir)

# Load and run the tests
loader = unittest.TestLoader()
tests = loader.discover(
    os.path.join(os.getcwd(), "test", "unit"), top_level_dir=os.getcwd()
)
testRunner = unittest.runner.TextTestRunner()
testRunner.run(tests)
