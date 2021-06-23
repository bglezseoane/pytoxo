#!/usr/bin/env python

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

"""Run all accuracy test helper script.

Simple script to discover and run test inside the current directory
with the project home directory as working directory."""

import os
import unittest

# Assert project home directory because tests have routes since here
if os.path.basename(os.getcwd()) == "test":
    os.chdir(os.pardir)

# Load and run the tests
loader = unittest.TestLoader()
tests = loader.discover(
    os.path.join(os.getcwd(), "test", "accuracy"), top_level_dir=os.getcwd()
)
testRunner = unittest.runner.TextTestRunner()
testRunner.run(tests)
