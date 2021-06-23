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

"""PyToxo model unit test suite."""

import os
import tempfile
import unittest

import numpy
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

    def test_dict_parsing(self):
        """Test model genotypes dictionary parsing."""
        m = pytoxo.model.Model(
            genotypes_dict={
                "AABBCC": "x",
                "AABBCc": "x*(1+y)",
                "AABBcc": "x*(1+y)^2",
                "AABbCC": "x*(1+y)",
                "AABbCc": "x*(1+y)^2",
                "AABbcc": "x*(1+y)^3",
                "AAbbCC": "x*(1+y)^2",
                "AAbbCc": "x*(1+y)^3",
                "AAbbcc": "x*(1+y)^4",
                "AaBBCC": "x*(1+y)",
                "AaBBCc": "x*(1+y)^2",
                "AaBBcc": "x*(1+y)^3",
                "AaBbCC": "x*(1+y)^2",
                "AaBbCc": "x*(1+y)^3",
                "AaBbcc": "x*(1+y)^4",
                "AabbCC": "x*(1+y)^3",
                "AabbCc": "x*(1+y)^4",
                "Aabbcc": "x*(1+y)^5",
                "aaBBCC": "x*(1+y)^2",
                "aaBBCc": "x*(1+y)^3",
                "aaBBcc": "x*(1+y)^4",
                "aaBbCC": "x*(1+y)^3",
                "aaBbCc": "x*(1+y)^4",
                "aaBbcc": "x*(1+y)^5",
                "aabbCC": "x*(1+y)^4",
                "aabbCc": "x*(1+y)^5",
                "aabbcc": "x*(1+y)^6",
            }  # Based in `additive_3` model
        )

        # Name
        self.assertEqual("unnamed", m._name)

        # Order
        self.assertEqual(3, m._order)

        # Penetrances
        """The following expressions are those of the file
        `models/additive_3.csv`, but corrected with small modifications
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
        self.assertEqual("x*(y + 1)", str(m._penetrances[9]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[10]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[11]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[12]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[13]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[14]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[15]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[16]))
        self.assertEqual("x*(y + 1)**5", str(m._penetrances[17]))
        self.assertEqual("x*(y + 1)**2", str(m._penetrances[18]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[19]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[20]))
        self.assertEqual("x*(y + 1)**3", str(m._penetrances[21]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[22]))
        self.assertEqual("x*(y + 1)**5", str(m._penetrances[23]))
        self.assertEqual("x*(y + 1)**4", str(m._penetrances[24]))
        self.assertEqual("x*(y + 1)**5", str(m._penetrances[25]))
        self.assertEqual("x*(y + 1)**6", str(m._penetrances[26]))
        # Some type checks
        self.assertEqual(sympy.Symbol, type(m._penetrances[0]))
        self.assertEqual(sympy.Mul, type(m._penetrances[5]))
        self.assertEqual(sympy.Mul, type(m._penetrances[6]))
        self.assertEqual(sympy.Mul, type(m._penetrances[8]))
        self.assertEqual(sympy.Mul, type(m._penetrances[12]))
        self.assertEqual(sympy.Mul, type(m._penetrances[26]))

        # Variables
        self.assertEqual([sympy.abc.x, sympy.abc.y], m._variables)

    def test_all_parsing_modes_equality(self):
        """Test that the three parsing modes produce equivalent model objects,
        given equivalent input configurations."""
        self.assertEqual(
            pytoxo.model.Model(
                filename=os.path.join("models", "additive_3.csv"), model_name="testing"
            ),
            pytoxo.model.Model(
                genotypes_dict={
                    "AABBCC": "x",
                    "AABBCc": "x*(1+y)",
                    "AABBcc": "x*(1+y)^2",
                    "AABbCC": "x*(1+y)",
                    "AABbCc": "x*(1+y)^2",
                    "AABbcc": "x*(1+y)^3",
                    "AAbbCC": "x*(1+y)^2",
                    "AAbbCc": "x*(1+y)^3",
                    "AAbbcc": "x*(1+y)^4",
                    "AaBBCC": "x*(1+y)",
                    "AaBBCc": "x*(1+y)^2",
                    "AaBBcc": "x*(1+y)^3",
                    "AaBbCC": "x*(1+y)^2",
                    "AaBbCc": "x*(1+y)^3",
                    "AaBbcc": "x*(1+y)^4",
                    "AabbCC": "x*(1+y)^3",
                    "AabbCc": "x*(1+y)^4",
                    "Aabbcc": "x*(1+y)^5",
                    "aaBBCC": "x*(1+y)^2",
                    "aaBBCc": "x*(1+y)^3",
                    "aaBBcc": "x*(1+y)^4",
                    "aaBbCC": "x*(1+y)^3",
                    "aaBbCc": "x*(1+y)^4",
                    "aaBbcc": "x*(1+y)^5",
                    "aabbCC": "x*(1+y)^4",
                    "aabbCc": "x*(1+y)^5",
                    "aabbcc": "x*(1+y)^6",
                },
                model_name="testing",
            ),
            pytoxo.model.Model(
                definitions=numpy.array(
                    [
                        "AABBCC",
                        "AABBCc",
                        "AABBcc",
                        "AABbCC",
                        "AABbCc",
                        "AABbcc",
                        "AAbbCC",
                        "AAbbCc",
                        "AAbbcc",
                        "AaBBCC",
                        "AaBBCc",
                        "AaBBcc",
                        "AaBbCC",
                        "AaBbCc",
                        "AaBbcc",
                        "AabbCC",
                        "AabbCc",
                        "Aabbcc",
                        "aaBBCC",
                        "aaBBCc",
                        "aaBBcc",
                        "aaBbCC",
                        "aaBbCc",
                        "aaBbcc",
                        "aabbCC",
                        "aabbCc",
                        "aabbcc",
                    ]
                ),
                probabilities=numpy.array(
                    [
                        "x",
                        "x*(1+y)",
                        "x*(1+y)^2",
                        "x*(1+y)",
                        "x*(1+y)^2",
                        "x*(1+y)^3",
                        "x*(1+y)^2",
                        "x*(1+y)^3",
                        "x*(1+y)^4",
                        "x*(1+y)",
                        "x*(1+y)^2",
                        "x*(1+y)^3",
                        "x*(1+y)^2",
                        "x*(1+y)^3",
                        "x*(1+y)^4",
                        "x*(1+y)^3",
                        "x*(1+y)^4",
                        "x*(1+y)^5",
                        "x*(1+y)^2",
                        "x*(1+y)^3",
                        "x*(1+y)^4",
                        "x*(1+y)^3",
                        "x*(1+y)^4",
                        "x*(1+y)^5",
                        "x*(1+y)^4",
                        "x*(1+y)^5",
                        "x*(1+y)^6",
                    ]
                ),
                model_name="testing",
            ),
        )

        self.assertEqual(
            pytoxo.model.Model(
                genotypes_dict={
                    "AABBCCDDEE": "x",
                    "AABBCCDDEe": "x",
                    "AABBCCDDee": "x",
                    "AABBCCDdEE": "x",
                    "AABBCCDdEe": "x",
                    "AABBCCDdee": "x",
                    "AABBCCddEE": "x",
                    "AABBCCddEe": "x",
                    "AABBCCddee": "x",
                    "AABBCcDDEE": "x",
                    "AABBCcDDEe": "x",
                    "AABBCcDDee": "x",
                    "AABBCcDdEE": "x",
                    "AABBCcDdEe": "x",
                    "AABBCcDdee": "x",
                    "AABBCcddEE": "x",
                    "AABBCcddEe": "x",
                    "AABBCcddee": "x",
                    "AABBccDDEE": "x",
                    "AABBccDDEe": "x",
                    "AABBccDDee": "x",
                    "AABBccDdEE": "x",
                    "AABBccDdEe": "x",
                    "AABBccDdee": "x",
                    "AABBccddEE": "x",
                    "AABBccddEe": "x",
                    "AABBccddee": "x",
                    "AABbCCDDEE": "x",
                    "AABbCCDDEe": "x",
                    "AABbCCDDee": "x",
                    "AABbCCDdEE": "x",
                    "AABbCCDdEe": "x",
                    "AABbCCDdee": "x",
                    "AABbCCddEE": "x",
                    "AABbCCddEe": "x",
                    "AABbCCddee": "x",
                    "AABbCcDDEE": "x",
                    "AABbCcDDEe": "x",
                    "AABbCcDDee": "x",
                    "AABbCcDdEE": "x",
                    "AABbCcDdEe": "x",
                    "AABbCcDdee": "x",
                    "AABbCcddEE": "x",
                    "AABbCcddEe": "x",
                    "AABbCcddee": "x",
                    "AABbccDDEE": "x",
                    "AABbccDDEe": "x",
                    "AABbccDDee": "x",
                    "AABbccDdEE": "x",
                    "AABbccDdEe": "x",
                    "AABbccDdee": "x",
                    "AABbccddEE": "x",
                    "AABbccddEe": "x",
                    "AABbccddee": "x",
                    "AAbbCCDDEE": "x",
                    "AAbbCCDDEe": "x",
                    "AAbbCCDDee": "x",
                    "AAbbCCDdEE": "x",
                    "AAbbCCDdEe": "x",
                    "AAbbCCDdee": "x",
                    "AAbbCCddEE": "x",
                    "AAbbCCddEe": "x",
                    "AAbbCCddee": "x",
                    "AAbbCcDDEE": "x",
                    "AAbbCcDDEe": "x",
                    "AAbbCcDDee": "x",
                    "AAbbCcDdEE": "x",
                    "AAbbCcDdEe": "x",
                    "AAbbCcDdee": "x",
                    "AAbbCcddEE": "x",
                    "AAbbCcddEe": "x",
                    "AAbbCcddee": "x",
                    "AAbbccDDEE": "x",
                    "AAbbccDDEe": "x",
                    "AAbbccDDee": "x",
                    "AAbbccDdEE": "x",
                    "AAbbccDdEe": "x",
                    "AAbbccDdee": "x",
                    "AAbbccddEE": "x",
                    "AAbbccddEe": "x",
                    "AAbbccddee": "x",
                    "AaBBCCDDEE": "x",
                    "AaBBCCDDEe": "x",
                    "AaBBCCDDee": "x",
                    "AaBBCCDdEE": "x",
                    "AaBBCCDdEe": "x",
                    "AaBBCCDdee": "x",
                    "AaBBCCddEE": "x",
                    "AaBBCCddEe": "x",
                    "AaBBCCddee": "x",
                    "AaBBCcDDEE": "x",
                    "AaBBCcDDEe": "x",
                    "AaBBCcDDee": "x",
                    "AaBBCcDdEE": "x",
                    "AaBBCcDdEe": "x",
                    "AaBBCcDdee": "x",
                    "AaBBCcddEE": "x",
                    "AaBBCcddEe": "x",
                    "AaBBCcddee": "x",
                    "AaBBccDDEE": "x",
                    "AaBBccDDEe": "x",
                    "AaBBccDDee": "x",
                    "AaBBccDdEE": "x",
                    "AaBBccDdEe": "x",
                    "AaBBccDdee": "x",
                    "AaBBccddEE": "x",
                    "AaBBccddEe": "x",
                    "AaBBccddee": "x",
                    "AaBbCCDDEE": "x",
                    "AaBbCCDDEe": "x",
                    "AaBbCCDDee": "x",
                    "AaBbCCDdEE": "x",
                    "AaBbCCDdEe": "x",
                    "AaBbCCDdee": "x",
                    "AaBbCCddEE": "x",
                    "AaBbCCddEe": "x",
                    "AaBbCCddee": "x",
                    "AaBbCcDDEE": "x",
                    "AaBbCcDDEe": "x",
                    "AaBbCcDDee": "x",
                    "AaBbCcDdEE": "x",
                    "AaBbCcDdEe": "x*(1+y)",
                    "AaBbCcDdee": "x*(1+y)^2",
                    "AaBbCcddEE": "x",
                    "AaBbCcddEe": "x*(1+y)^2",
                    "AaBbCcddee": "x*(1+y)^4",
                    "AaBbccDDEE": "x",
                    "AaBbccDDEe": "x",
                    "AaBbccDDee": "x",
                    "AaBbccDdEE": "x",
                    "AaBbccDdEe": "x*(1+y)^2",
                    "AaBbccDdee": "x*(1+y)^4",
                    "AaBbccddEE": "x",
                    "AaBbccddEe": "x*(1+y)^4",
                    "AaBbccddee": "x*(1+y)^8",
                    "AabbCCDDEE": "x",
                    "AabbCCDDEe": "x",
                    "AabbCCDDee": "x",
                    "AabbCCDdEE": "x",
                    "AabbCCDdEe": "x",
                    "AabbCCDdee": "x",
                    "AabbCCddEE": "x",
                    "AabbCCddEe": "x",
                    "AabbCCddee": "x",
                    "AabbCcDDEE": "x",
                    "AabbCcDDEe": "x",
                    "AabbCcDDee": "x",
                    "AabbCcDdEE": "x",
                    "AabbCcDdEe": "x*(1+y)^2",
                    "AabbCcDdee": "x*(1+y)^4",
                    "AabbCcddEE": "x",
                    "AabbCcddEe": "x*(1+y)^4",
                    "AabbCcddee": "x*(1+y)^8",
                    "AabbccDDEE": "x",
                    "AabbccDDEe": "x",
                    "AabbccDDee": "x",
                    "AabbccDdEE": "x",
                    "AabbccDdEe": "x*(1+y)^4",
                    "AabbccDdee": "x*(1+y)^8",
                    "AabbccddEE": "x",
                    "AabbccddEe": "x*(1+y)^8",
                    "Aabbccddee": "x*(1+y)^16",
                    "aaBBCCDDEE": "x",
                    "aaBBCCDDEe": "x",
                    "aaBBCCDDee": "x",
                    "aaBBCCDdEE": "x",
                    "aaBBCCDdEe": "x",
                    "aaBBCCDdee": "x",
                    "aaBBCCddEE": "x",
                    "aaBBCCddEe": "x",
                    "aaBBCCddee": "x",
                    "aaBBCcDDEE": "x",
                    "aaBBCcDDEe": "x",
                    "aaBBCcDDee": "x",
                    "aaBBCcDdEE": "x",
                    "aaBBCcDdEe": "x",
                    "aaBBCcDdee": "x",
                    "aaBBCcddEE": "x",
                    "aaBBCcddEe": "x",
                    "aaBBCcddee": "x",
                    "aaBBccDDEE": "x",
                    "aaBBccDDEe": "x",
                    "aaBBccDDee": "x",
                    "aaBBccDdEE": "x",
                    "aaBBccDdEe": "x",
                    "aaBBccDdee": "x",
                    "aaBBccddEE": "x",
                    "aaBBccddEe": "x",
                    "aaBBccddee": "x",
                    "aaBbCCDDEE": "x",
                    "aaBbCCDDEe": "x",
                    "aaBbCCDDee": "x",
                    "aaBbCCDdEE": "x",
                    "aaBbCCDdEe": "x",
                    "aaBbCCDdee": "x",
                    "aaBbCCddEE": "x",
                    "aaBbCCddEe": "x",
                    "aaBbCCddee": "x",
                    "aaBbCcDDEE": "x",
                    "aaBbCcDDEe": "x",
                    "aaBbCcDDee": "x",
                    "aaBbCcDdEE": "x",
                    "aaBbCcDdEe": "x*(1+y)^2",
                    "aaBbCcDdee": "x*(1+y)^4",
                    "aaBbCcddEE": "x",
                    "aaBbCcddEe": "x*(1+y)^4",
                    "aaBbCcddee": "x*(1+y)^8",
                    "aaBbccDDEE": "x",
                    "aaBbccDDEe": "x",
                    "aaBbccDDee": "x",
                    "aaBbccDdEE": "x",
                    "aaBbccDdEe": "x*(1+y)^4",
                    "aaBbccDdee": "x*(1+y)^8",
                    "aaBbccddEE": "x",
                    "aaBbccddEe": "x*(1+y)^8",
                    "aaBbccddee": "x*(1+y)^16",
                    "aabbCCDDEE": "x",
                    "aabbCCDDEe": "x",
                    "aabbCCDDee": "x",
                    "aabbCCDdEE": "x",
                    "aabbCCDdEe": "x",
                    "aabbCCDdee": "x",
                    "aabbCCddEE": "x",
                    "aabbCCddEe": "x",
                    "aabbCCddee": "x",
                    "aabbCcDDEE": "x",
                    "aabbCcDDEe": "x",
                    "aabbCcDDee": "x",
                    "aabbCcDdEE": "x",
                    "aabbCcDdEe": "x*(1+y)^4",
                    "aabbCcDdee": "x*(1+y)^8",
                    "aabbCcddEE": "x",
                    "aabbCcddEe": "x*(1+y)^8",
                    "aabbCcddee": "x*(1+y)^16",
                    "aabbccDDEE": "x",
                    "aabbccDDEe": "x",
                    "aabbccDDee": "x",
                    "aabbccDdEE": "x",
                    "aabbccDdEe": "x*(1+y)^8",
                    "aabbccDdee": "x*(1+y)^16",
                    "aabbccddEE": "x",
                    "aabbccddEe": "x*(1+y)^16",
                    "aabbccddee": "x*(1+y)^32",
                },
                model_name="testing2",
            ),
            pytoxo.model.Model(
                filename=os.path.join("models", "multiplicative_5.csv"),
                model_name="testing2",
            ),
            pytoxo.model.Model(
                definitions=[
                    "AABBCCDDEE",
                    "AABBCCDDEe",
                    "AABBCCDDee",
                    "AABBCCDdEE",
                    "AABBCCDdEe",
                    "AABBCCDdee",
                    "AABBCCddEE",
                    "AABBCCddEe",
                    "AABBCCddee",
                    "AABBCcDDEE",
                    "AABBCcDDEe",
                    "AABBCcDDee",
                    "AABBCcDdEE",
                    "AABBCcDdEe",
                    "AABBCcDdee",
                    "AABBCcddEE",
                    "AABBCcddEe",
                    "AABBCcddee",
                    "AABBccDDEE",
                    "AABBccDDEe",
                    "AABBccDDee",
                    "AABBccDdEE",
                    "AABBccDdEe",
                    "AABBccDdee",
                    "AABBccddEE",
                    "AABBccddEe",
                    "AABBccddee",
                    "AABbCCDDEE",
                    "AABbCCDDEe",
                    "AABbCCDDee",
                    "AABbCCDdEE",
                    "AABbCCDdEe",
                    "AABbCCDdee",
                    "AABbCCddEE",
                    "AABbCCddEe",
                    "AABbCCddee",
                    "AABbCcDDEE",
                    "AABbCcDDEe",
                    "AABbCcDDee",
                    "AABbCcDdEE",
                    "AABbCcDdEe",
                    "AABbCcDdee",
                    "AABbCcddEE",
                    "AABbCcddEe",
                    "AABbCcddee",
                    "AABbccDDEE",
                    "AABbccDDEe",
                    "AABbccDDee",
                    "AABbccDdEE",
                    "AABbccDdEe",
                    "AABbccDdee",
                    "AABbccddEE",
                    "AABbccddEe",
                    "AABbccddee",
                    "AAbbCCDDEE",
                    "AAbbCCDDEe",
                    "AAbbCCDDee",
                    "AAbbCCDdEE",
                    "AAbbCCDdEe",
                    "AAbbCCDdee",
                    "AAbbCCddEE",
                    "AAbbCCddEe",
                    "AAbbCCddee",
                    "AAbbCcDDEE",
                    "AAbbCcDDEe",
                    "AAbbCcDDee",
                    "AAbbCcDdEE",
                    "AAbbCcDdEe",
                    "AAbbCcDdee",
                    "AAbbCcddEE",
                    "AAbbCcddEe",
                    "AAbbCcddee",
                    "AAbbccDDEE",
                    "AAbbccDDEe",
                    "AAbbccDDee",
                    "AAbbccDdEE",
                    "AAbbccDdEe",
                    "AAbbccDdee",
                    "AAbbccddEE",
                    "AAbbccddEe",
                    "AAbbccddee",
                    "AaBBCCDDEE",
                    "AaBBCCDDEe",
                    "AaBBCCDDee",
                    "AaBBCCDdEE",
                    "AaBBCCDdEe",
                    "AaBBCCDdee",
                    "AaBBCCddEE",
                    "AaBBCCddEe",
                    "AaBBCCddee",
                    "AaBBCcDDEE",
                    "AaBBCcDDEe",
                    "AaBBCcDDee",
                    "AaBBCcDdEE",
                    "AaBBCcDdEe",
                    "AaBBCcDdee",
                    "AaBBCcddEE",
                    "AaBBCcddEe",
                    "AaBBCcddee",
                    "AaBBccDDEE",
                    "AaBBccDDEe",
                    "AaBBccDDee",
                    "AaBBccDdEE",
                    "AaBBccDdEe",
                    "AaBBccDdee",
                    "AaBBccddEE",
                    "AaBBccddEe",
                    "AaBBccddee",
                    "AaBbCCDDEE",
                    "AaBbCCDDEe",
                    "AaBbCCDDee",
                    "AaBbCCDdEE",
                    "AaBbCCDdEe",
                    "AaBbCCDdee",
                    "AaBbCCddEE",
                    "AaBbCCddEe",
                    "AaBbCCddee",
                    "AaBbCcDDEE",
                    "AaBbCcDDEe",
                    "AaBbCcDDee",
                    "AaBbCcDdEE",
                    "AaBbCcDdEe",
                    "AaBbCcDdee",
                    "AaBbCcddEE",
                    "AaBbCcddEe",
                    "AaBbCcddee",
                    "AaBbccDDEE",
                    "AaBbccDDEe",
                    "AaBbccDDee",
                    "AaBbccDdEE",
                    "AaBbccDdEe",
                    "AaBbccDdee",
                    "AaBbccddEE",
                    "AaBbccddEe",
                    "AaBbccddee",
                    "AabbCCDDEE",
                    "AabbCCDDEe",
                    "AabbCCDDee",
                    "AabbCCDdEE",
                    "AabbCCDdEe",
                    "AabbCCDdee",
                    "AabbCCddEE",
                    "AabbCCddEe",
                    "AabbCCddee",
                    "AabbCcDDEE",
                    "AabbCcDDEe",
                    "AabbCcDDee",
                    "AabbCcDdEE",
                    "AabbCcDdEe",
                    "AabbCcDdee",
                    "AabbCcddEE",
                    "AabbCcddEe",
                    "AabbCcddee",
                    "AabbccDDEE",
                    "AabbccDDEe",
                    "AabbccDDee",
                    "AabbccDdEE",
                    "AabbccDdEe",
                    "AabbccDdee",
                    "AabbccddEE",
                    "AabbccddEe",
                    "Aabbccddee",
                    "aaBBCCDDEE",
                    "aaBBCCDDEe",
                    "aaBBCCDDee",
                    "aaBBCCDdEE",
                    "aaBBCCDdEe",
                    "aaBBCCDdee",
                    "aaBBCCddEE",
                    "aaBBCCddEe",
                    "aaBBCCddee",
                    "aaBBCcDDEE",
                    "aaBBCcDDEe",
                    "aaBBCcDDee",
                    "aaBBCcDdEE",
                    "aaBBCcDdEe",
                    "aaBBCcDdee",
                    "aaBBCcddEE",
                    "aaBBCcddEe",
                    "aaBBCcddee",
                    "aaBBccDDEE",
                    "aaBBccDDEe",
                    "aaBBccDDee",
                    "aaBBccDdEE",
                    "aaBBccDdEe",
                    "aaBBccDdee",
                    "aaBBccddEE",
                    "aaBBccddEe",
                    "aaBBccddee",
                    "aaBbCCDDEE",
                    "aaBbCCDDEe",
                    "aaBbCCDDee",
                    "aaBbCCDdEE",
                    "aaBbCCDdEe",
                    "aaBbCCDdee",
                    "aaBbCCddEE",
                    "aaBbCCddEe",
                    "aaBbCCddee",
                    "aaBbCcDDEE",
                    "aaBbCcDDEe",
                    "aaBbCcDDee",
                    "aaBbCcDdEE",
                    "aaBbCcDdEe",
                    "aaBbCcDdee",
                    "aaBbCcddEE",
                    "aaBbCcddEe",
                    "aaBbCcddee",
                    "aaBbccDDEE",
                    "aaBbccDDEe",
                    "aaBbccDDee",
                    "aaBbccDdEE",
                    "aaBbccDdEe",
                    "aaBbccDdee",
                    "aaBbccddEE",
                    "aaBbccddEe",
                    "aaBbccddee",
                    "aabbCCDDEE",
                    "aabbCCDDEe",
                    "aabbCCDDee",
                    "aabbCCDdEE",
                    "aabbCCDdEe",
                    "aabbCCDdee",
                    "aabbCCddEE",
                    "aabbCCddEe",
                    "aabbCCddee",
                    "aabbCcDDEE",
                    "aabbCcDDEe",
                    "aabbCcDDee",
                    "aabbCcDdEE",
                    "aabbCcDdEe",
                    "aabbCcDdee",
                    "aabbCcddEE",
                    "aabbCcddEe",
                    "aabbCcddee",
                    "aabbccDDEE",
                    "aabbccDDEe",
                    "aabbccDDee",
                    "aabbccDdEE",
                    "aabbccDdEe",
                    "aabbccDdee",
                    "aabbccddEE",
                    "aabbccddEe",
                    "aabbccddee",
                ],
                probabilities=[
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)",
                    "x*(1+y)^2",
                    "x",
                    "x*(1+y)^2",
                    "x*(1+y)^4",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^2",
                    "x*(1+y)^4",
                    "x",
                    "x*(1+y)^4",
                    "x*(1+y)^8",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^2",
                    "x*(1+y)^4",
                    "x",
                    "x*(1+y)^4",
                    "x*(1+y)^8",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^4",
                    "x*(1+y)^8",
                    "x",
                    "x*(1+y)^8",
                    "x*(1+y)^16",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^2",
                    "x*(1+y)^4",
                    "x",
                    "x*(1+y)^4",
                    "x*(1+y)^8",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^4",
                    "x*(1+y)^8",
                    "x",
                    "x*(1+y)^8",
                    "x*(1+y)^16",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^4",
                    "x*(1+y)^8",
                    "x",
                    "x*(1+y)^8",
                    "x*(1+y)^16",
                    "x",
                    "x",
                    "x",
                    "x",
                    "x*(1+y)^8",
                    "x*(1+y)^16",
                    "x",
                    "x*(1+y)^16",
                    "x*(1+y)^32",
                ],
                model_name="testing2",
            ),
        )

    def test_file_parsing_different_var_names(self):
        """Test model files parsing. This version uses unusual names for the
        variables and not typical `x` and `y` to assert function."""
        with tempfile.TemporaryDirectory() as tmp_root:
            # Use a real model as base to the test
            with open(os.path.join("models", "additive_2.csv"), "r") as f:
                additive_2_model_content = f.read()
            # Substitute original model `x` and `y` with `g` and `w`, respectively
            additive_2_model_content = additive_2_model_content.replace("x", "g")
            additive_2_model_content = additive_2_model_content.replace("y", "w")
            # Save to a temp file
            tmp_file = os.path.join(tmp_root, "additive_2.csv")
            with open(tmp_file, "x") as f:
                f.write(additive_2_model_content)

            # Read the modified model and go on with the test...
            m = pytoxo.model.Model(tmp_file)

            # Name
            self.assertEqual("additive_2", m._name)

            # Order
            self.assertEqual(2, m._order)

            # Penetrances
            """The following expressions are those of the file 
            `models/additive_2.csv`, but corrected with small modifications
            that the library Sympy does when normalizing them. And with 
            the modifications `x` and `y` to `g` and `w`, respectively"""
            self.assertEqual("g", str(m._penetrances[0]))
            self.assertEqual("g*(w + 1)", str(m._penetrances[1]))
            self.assertEqual("g*(w + 1)**2", str(m._penetrances[2]))
            self.assertEqual("g*(w + 1)", str(m._penetrances[3]))
            self.assertEqual("g*(w + 1)**2", str(m._penetrances[4]))
            self.assertEqual("g*(w + 1)**3", str(m._penetrances[5]))
            self.assertEqual("g*(w + 1)**2", str(m._penetrances[6]))
            self.assertEqual("g*(w + 1)**3", str(m._penetrances[7]))
            self.assertEqual("g*(w + 1)**4", str(m._penetrances[8]))
            # Some type checks
            self.assertEqual(sympy.Symbol, type(m._penetrances[0]))
            self.assertEqual(sympy.Mul, type(m._penetrances[1]))
            self.assertEqual(sympy.Mul, type(m._penetrances[2]))
            self.assertEqual(sympy.Mul, type(m._penetrances[3]))
            self.assertEqual(sympy.Mul, type(m._penetrances[4]))
            self.assertEqual(sympy.Mul, type(m._penetrances[5]))

            # Variables
            self.assertEqual([sympy.abc.g, sympy.abc.w], m._variables)

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
                    pytoxo.errors.BadFormedModelError,
                    lambda: pytoxo.model.Model(bad_formed_models_path),
                )

    def test_bad_init_param_use_raising(self):
        """Test raising during `Model` init, due to a bad use of the
        parameters."""
        with self.assertRaises(ValueError):
            _ = pytoxo.Model(definitions=["AABb", "AABB"])
        with self.assertRaises(ValueError):
            _ = pytoxo.Model(probabilities=["x", "y"])
        with self.assertRaises(ValueError):
            _ = pytoxo.Model(
                filename=os.path.join("models", "multiplicative_4.csv"),
                model_name="error",
                definitions=["AABb", "AABB"],
            )
        with self.assertRaises(ValueError):
            _ = pytoxo.Model(
                filename=os.path.join("models", "multiplicative_4.csv"),
                model_name="error",
                probabilities=["x", "y"],
            )
        with self.assertRaises(ValueError):
            _ = pytoxo.Model(
                filename=os.path.join("models", "multiplicative_4.csv"),
                model_name="error",
                probabilities=["x", "y"],
                definitions=["AABb", "AABB"],
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

    def test_check_solution(self):
        """Test model `_check_solution` method."""
        m = pytoxo.model.Model(
            os.path.join("models", "multiplicative_2.csv")
        )  # The only relevant detail to this test is that variables are `x` and `y`

        eqs1 = [
            sympy.Eq(sympy.abc.x + 2, 3),
            sympy.Eq(sympy.abc.x ** 2 + sympy.abc.x ** 2, 2),
        ]
        self.assertTrue(
            m._check_solution(eqs1, ({sympy.abc.x: 1.0, sympy.abc.y: 0}))[0]
        )

        eqs2 = [
            sympy.Eq(sympy.abc.x ** 2 + sympy.abc.y, 4),
            sympy.Eq(sympy.abc.y, 0),
        ]
        self.assertTrue(
            m._check_solution(eqs2, ({sympy.abc.x: 2.0, sympy.abc.y: 0}))[0]
        )
        self.assertFalse(
            m._check_solution(eqs2, ({sympy.abc.x: 2.0, sympy.abc.y: 1}))[0]
        )
        self.assertTrue(
            m._check_solution(
                eqs2,
                (
                    {
                        sympy.abc.x: 2.0,
                        sympy.abc.y: m.calculate_tolerable_solution_error_delta(),
                    }
                ),
            )[0]
        )
        self.assertFalse(
            m._check_solution(
                eqs2,
                (
                    {
                        sympy.abc.x: 2.0,
                        sympy.abc.y: m.calculate_tolerable_solution_error_delta() * 10,
                    }
                ),
            )[0]
        )
        self.assertFalse(
            m._check_solution(
                eqs2,
                (
                    {
                        sympy.abc.x: 2.0 + m.calculate_tolerable_solution_error_delta(),
                        sympy.abc.y: 0,
                    }
                ),
            )[0]
        )

    def test_model_comparison(self):
        """Equality criteria between model objects."""
        m1 = pytoxo.model.Model(filename=os.path.join("models", "multiplicative_2.csv"))
        m2 = pytoxo.model.Model(filename=os.path.join("models", "multiplicative_2.csv"))
        m3 = pytoxo.model.Model(filename=os.path.join("models", "multiplicative_3.csv"))
        m4 = pytoxo.model.Model(
            filename=os.path.join("models", "multiplicative_3.csv"),
            model_name="other_name",
        )
        self.assertEqual(m1, m1)
        self.assertEqual(m1, m2)
        self.assertNotEqual(m2, m3)
        self.assertNotEqual(m1, m3)
        self.assertNotEqual(m3, m4)
        self.assertEqual(m4, m4)

    def test_probabilities_sort(self):
        """Test that the probability expressions are sorted attending to
        genotype definitions alphabetical sort. This is necessary to assert the association between
        genotype definitions and probabilities during the calculus process
        and in the final penetrance table."""
        m1 = pytoxo.model.Model(
            genotypes_dict={
                "AaBbcc": "x*(1+y)^4",
                "AabbCC": "x*(1+y)^3",
                "AabbCc": "x*(1+y)^4",
                "Aabbcc": "x*(1+y)^5",
                "aaBBCC": "x*(1+y)^2",
                "aaBBCc": "x*(1+y)^3",
                "aaBBcc": "x*(1+y)^4",
                "AABbCC": "x*(1+y)",
                "AABbCc": "x*(1+y)^2",
                "AABbcc": "x*(1+y)^3",
                "AAbbCC": "x*(1+y)^2",
                "AAbbCc": "x*(1+y)^3",
                "AAbbcc": "x*(1+y)^4",
                "aaBbCC": "x*(1+y)^3",
                "aaBbCc": "x*(1+y)^4",
                "aaBbcc": "x*(1+y)^5",
                "aabbCC": "x*(1+y)^4",
                "aabbCc": "x*(1+y)^5",
                "aabbcc": "x*(1+y)^6",
                "AABBCC": "x",
                "AABBCc": "x*(1+y)",
                "AABBcc": "x*(1+y)^2",
                "AaBBCC": "x*(1+y)",
                "AaBBCc": "x*(1+y)^2",
                "AaBBcc": "x*(1+y)^3",
                "AaBbCC": "x*(1+y)^2",
                "AaBbCc": "x*(1+y)^3",
            },
            model_name="additive_3",  # Based in `additive_3` model, unsorted
        )

        m2 = pytoxo.model.Model(
            filename=os.path.join("models", "additive_3.csv")
        )  # This one is already sorted

        self.assertEqual(m1, m2)

        # Also check some probability expressions manually
        self.assertEqual("x", str(m1._penetrances[0]))
        self.assertEqual("x*(y + 1)", str(m1._penetrances[1]))
        self.assertEqual("x*(y + 1)**2", str(m2._penetrances[2]))
        self.assertEqual("x*(y + 1)", str(m2._penetrances[3]))
        self.assertEqual("x*(y + 1)**2", str(m2._penetrances[4]))
        self.assertEqual("x*(y + 1)**3", str(m2._penetrances[5]))
        self.assertEqual("x*(y + 1)**2", str(m2._penetrances[12]))
        self.assertEqual("x*(y + 1)**3", str(m2._penetrances[13]))
        self.assertEqual("x*(y + 1)**4", str(m2._penetrances[14]))
        self.assertEqual("x*(y + 1)**3", str(m2._penetrances[15]))
        self.assertEqual("x*(y + 1)**3", str(m2._penetrances[19]))
        self.assertEqual("x*(y + 1)**4", str(m2._penetrances[20]))
        self.assertEqual("x*(y + 1)**3", str(m2._penetrances[21]))
        self.assertEqual("x*(y + 1)**6", str(m2._penetrances[26]))

    def test_find_parameters_check(self):
        m = pytoxo.model.Model(filename=os.path.join("models", "multiplicative_2.csv"))

        # Only one MAF, two needed
        with self.assertRaises(ValueError):
            m.find_max_heritability_table([0.2], 0.4)

        # MAFs should be floats
        with self.assertRaises(ValueError):
            m.find_max_heritability_table([True, 0.3], 0.4)

        # MAFs should be less or equal than 0.5
        with self.assertRaises(ValueError):
            m.find_max_heritability_table([0.51] * m.order, 0.4)

        # Prevalence should be a float
        with self.assertRaises(ValueError):
            m.find_max_heritability_table([0.2] * m.order, 4)

        # Timeout should be an integer or a bool
        with self.assertRaises(ValueError):
            m.find_max_heritability_table([0.2] * m.order, 0.4, 1.5)

        # Check flag should be a bool
        with self.assertRaises(ValueError):
            m.find_max_heritability_table([0.2] * m.order, 0.4, 4, 1)
