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

"""Part of PyToxo integration test suite."""

import os
import unittest

import pytoxo.calculations
import pytoxo.errors
import pytoxo.model


class ErrorHandlingTestSuite(unittest.TestCase):
    """Test suite which checks the correct raising and catching of PyToxo
    derived exceptions.

    Tests for PyToxo full process at integration level.
    """

    def test_exception_raising_during_models_solving(self):
        """Test exception raising and catching during models solving using some
        known unsolvable cases."""
        # Cases known as unsolvable
        cases = [
            (os.path.join("models", "additive_4.csv"), 0.5, 0.0),
            (os.path.join("models", "additive_5.csv"), 0.4, 1.0),
            (os.path.join("models", "additive_6.csv"), 0.3, 1.0),
            (os.path.join("models", "additive_6.csv"), 0.2, 1.0),
        ]

        for case_file, case_mafs, case_p in cases:
            model = pytoxo.model.Model(case_file)
            with self.assertRaises(pytoxo.errors.UnsolvableModelError):
                model.find_max_heritability_table([case_mafs] * model.order, case_p)
