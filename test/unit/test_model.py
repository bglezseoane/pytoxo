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

"""PyToxo model unit test suite."""

import os
import tempfile
import unittest

import sympy
import sympy.abc

import pytoxo.errors
import pytoxo.model


class ModelUnitTestSuite(unittest.TestCase):
    """Tests for `pytoxo/model.py` at unit level.

    Warnings
    --------
    Tests are sensible to model files `models/additive_2.csv` and
    `models/multiplicative_2.csv`.
    """

    def test_file_parsing_1(self):
        """Test model files parsing."""
        m = pytoxo.model.Model(os.path.join("models", "multiplicative_2.csv"))

        # Name
        self.assertEqual("multiplicative_2", m._name)

        # Order
        self.assertEqual(2, m._order)

        # Penetrances
        """The following expressions are those of the file 
        `models/multiplicative_2.csv`, but corrected with small modifications
        that the library Sympy does when normalizing them. E.g .: symbol `**` 
        instead of `^`, spaces around additions, etc."""
        self.assertEqual("x", str(m._penetrances[0]))
        self.assertEqual("x", str(m._penetrances[1]))
        self.assertEqual("x", str(m._penetrances[2]))
        self.assertEqual("x", str(m._penetrances[3]))
        self.assertEqual("x*(y + 1)", str(m._penetrances[4]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[5]))
        self.assertEqual("x", str(m._penetrances[6]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[7]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[8]))
        # Some type checks
        self.assertEqual(sympy.Symbol, type(m._penetrances[0]))
        self.assertEqual(sympy.Mul, type(m._penetrances[4]))
        self.assertEqual(sympy.Mul, type(m._penetrances[8]))

        # Variables
        self.assertEqual([sympy.abc.x, sympy.abc.y], m._variables)

    def test_file_parsing_2(self):
        """Test model files parsing."""
        m = pytoxo.model.Model(os.path.join("models", "additive_2.csv"))

        # Name
        self.assertEqual("additive_2", m._name)

        # Order
        self.assertEqual(2, m._order)

        # Penetrances
        """The following expressions are those of the file 
        `models/additive_2.csv`, but corrected with small modifications
        that the library Sympy does when normalizing them. E.g .: symbol `**` 
        instead of `^`, spaces around additions, etc."""
        self.assertEqual("x", str(m._penetrances[0]))
        self.assertEqual("x*(y + 1)", str(m._penetrances[1]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[2]))
        self.assertEqual("x*(y + 1)", str(m._penetrances[3]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[4]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[5]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[6]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[7]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[8]))
        # Some type checks
        self.assertEqual(sympy.Symbol, type(m._penetrances[0]))
        self.assertEqual(sympy.Mul, type(m._penetrances[1]))
        self.assertEqual(sympy.Mul, type(m._penetrances[2]))
        self.assertEqual(sympy.Mul, type(m._penetrances[3]))
        self.assertEqual(sympy.Mul, type(m._penetrances[4]))
        self.assertEqual(sympy.Mul, type(m._penetrances[5]))

        # Variables
        self.assertEqual([sympy.abc.x, sympy.abc.y], m._variables)

    def test_file_parsing_error_detection(self):
        """Test error handling during model files parsing."""

        # Test nonexistent file raise
        self.assertRaises(OSError, lambda: pytoxo.model.Model("nonexistent_file.csv"))

        # Read a well formed sample model to corrupt it and assert error raises
        with open(os.path.join("models", "additive_2.csv"), "r") as f:
            well_formed_file_content = f.readlines()

        with tempfile.TemporaryDirectory() as mock_bad_formed_models_dir:
            # Create some bad formed models using well formed content
            bad_formed_models_paths = [
                os.path.join(mock_bad_formed_models_dir, filename)
                for filename in ["1.csv", "2.csv", "3.csv"]
            ]
            for bad_formed_models_path in bad_formed_models_paths:
                with open(bad_formed_models_path, "w+") as f:
                    for line in well_formed_file_content:
                        # Introduce some errors
                        if os.path.basename(bad_formed_models_path) == "1.csv":
                            if "AABB" in line:
                                line = line.replace("AABB", "AAB")
                        if os.path.basename(bad_formed_models_path) == "2.csv":
                            if line.startswith("A") or line.startswith("a"):
                                splitted_line = line.split(",")
                                line = f"{splitted_line[0][:3]}, {splitted_line[1]}"
                        if os.path.basename(bad_formed_models_path) == "3.csv":
                            if "x*(1+y)^3" in line:
                                line = line.replace("x*(1+y)^3", "x*+*(1+y)^3")
                        f.write(line)

            # Test bad formed files raise
            for bad_formed_models_path in bad_formed_models_paths:
                self.assertRaises(
                    pytoxo.errors.ModelCSVParseError,
                    lambda: pytoxo.model.Model(bad_formed_models_path),
                )

    def test_max_penetrance_1(self):
        """Test model `max_penetrance` method."""
        m = pytoxo.model.Model(os.path.join("models", "additive_3.csv"))

        """Output from Toxo's `m = max_penetrance(obj)` private method for 
        the given model, adapted to a Sympy object"""
        expected_output = sympy.sympify("x*(y + 1)**6")

        output = m._max_penetrance()

        self.assertEqual(expected_output, output)

    def test_max_penetrance_2(self):
        """Test model `max_penetrance` method."""
        m = pytoxo.model.Model(os.path.join("models", "multiplicative_4.csv"))

        """Output from Toxo's `m = max_penetrance(obj)` private method for 
        the given model, adapted to a Sympy object"""
        expected_output = sympy.sympify("x*(y + 1)**16")

        output = m._max_penetrance()

        self.assertEqual(expected_output, output)

    def test_max_penetrance_3(self):
        """Test model `max_penetrance` method.

        This test uses mock file content with the peculiarity that the
        polynomial that is known larger is not in last place, as it usually
        is in real models. This circumstance is verified because the way in
        which the loop of the `max_penetrance` function works.
        """
        mock_model_content = (
            "aabbCcddEeffGgHH, x\n"
            + "aabbCcddEeffGgHh, x*(1+y)^16\n"
            + "aabbCcddEeffGghh, x*(1+y)^32\n"
            + "aabbCcddEeffgghh, x*(1+y)^64\n"
            + "aabbCcddEeffggHH, x\n"
            + "aabbCcddEeffggHh, x*(1+y)^32\n"
        )

        with tempfile.TemporaryDirectory() as mock_model_dir:
            mock_model_file = os.path.join(mock_model_dir, "model.csv")
            with open(mock_model_file, "x") as f:
                f.write(mock_model_content)

            m = pytoxo.model.Model(mock_model_file)

            # Known larger polynomial as Sympy string
            expected_output = sympy.sympify("x*(y + 1)**64")

            output = m._max_penetrance()

            self.assertEqual(expected_output, output)
