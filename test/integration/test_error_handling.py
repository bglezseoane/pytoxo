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
        # Cases determined unsolvable by Toxo
        cases = [
            (os.path.join("models", "additive_4.csv"), 0.3, 0.1),
            (os.path.join("models", "additive_5.csv"), 0.3, 0.1),
            (os.path.join("models", "additive_6.csv"), 0.49, 0.1),
            (os.path.join("models", "additive_6.csv"), 0.49, 0.2),
            (os.path.join("models", "multiplicative_2.csv"), 0.1, 0.1),
            (os.path.join("models", "multiplicative_2.csv"), 0.2, 0.1),
            (os.path.join("models", "multiplicative_2.csv"), 0.3456, 0.2),
        ]

        for case_file, case_mafs, case_p in cases:
            model = pytoxo.model.Model(case_file)
            with self.assertRaises(pytoxo.errors.UnsolvableModelError):
                model.find_max_heritability_table([case_mafs] * model.order, case_p)
