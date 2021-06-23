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

"""Part of PyToxo integration test suite."""

import os
import sys
import unittest

import sympy

import pytoxo.calculations
import pytoxo.model


class EquationSubstitutionsTestSuite(unittest.TestCase):
    """Test suite which checks the calculation of the values of the variables
    used in a model after its solution. To do this, instead of using the Toxo
    contrast, like in other suites of PyToxo, replace the values of the
    variables in the initial equation again.

    Tests for PyToxo penetrance table calculations at integration level.
    """

    def test_equation_substitutions_additive_2(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_2"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_multiplicative_2(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_2"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_2(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_2"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_additive_3(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_3"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_multiplicative_3(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_3"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_3(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_3"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_additive_4(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_4"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_multiplicative_4(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_4"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_4(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_4"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_additive_5(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_5"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    @unittest.skip
    def test_equation_substitutions_multiplicative_5(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_5"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_5(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_5"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_additive_6(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_6"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    @unittest.skip
    def test_equation_substitutions_multiplicative_6(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_6"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_6(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_6"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_additive_7(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_7"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    @unittest.skip
    def test_equation_substitutions_multiplicative_7(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_7"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_7(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_7"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_additive_8(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "additive_8"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    @unittest.skip
    def test_equation_substitutions_multiplicative_8(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "multiplicative_8"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    def test_equation_substitutions_threshold_8(self):
        """Test the calculation of the values of the variables used in a model
        after its solution. To do this, replace the values of the variables
        in the initial equation again.

        The test check both heritability and prevalence maximizations.

        This test is for the `additive_2` model. MAF, heritability and
        prevalence values are equal for all tests.
        """
        model = "threshold_8"
        order = int(model.split("_")[1])

        # Run within the helper function
        self._helper_test_equation_substitutions(
            test=self,
            model_file=os.path.join("models", f"{model}.csv"),
            order=order,
            accuracy_places=sys.float_info.dig - order,
        )

    @staticmethod
    def _helper_test_equation_substitutions(
        test,
        model_file,
        order,
        accuracy_places,
    ):
        """Helper method with the test skeleton to the calculation of the
        solution of a model solutions of the `EquationSubstitutionsTestSuite`
        test suite.
        """
        mafs = []
        for _ in range(order):
            mafs.append(0.1)
        prevalence = 0.1
        heritability = 0.1

        # Generate the model
        model = pytoxo.model.Model(model_file)

        # Run the test maximizing prevalence and then heritability
        for her_or_prev, max_method in zip(
            [prevalence, heritability],
            [
                pytoxo.model.Model._build_max_prevalence_system,
                pytoxo.model.Model._build_max_heritability_system,
            ],
        ):
            eq_system = max_method(model, mafs, her_or_prev)
            eq1_lhs = eq_system[0].lhs  # Left hand side
            eq2_lhs = eq_system[1].lhs  # Left hand side

            # Calculate the solutions for the two variables
            vars_sol = model._solve(eq_system=eq_system)

            # Now substitute them in the system equations
            substituted_eq1_lhs = eq1_lhs.subs(
                {
                    model.variables[0]: vars_sol[model.variables[0]],
                    model.variables[1]: vars_sol[model.variables[1]],
                }
            )
            substituted_eq2_lhs = eq2_lhs.subs(
                {
                    model.variables[0]: vars_sol[model.variables[0]],
                    model.variables[1]: vars_sol[model.variables[1]],
                }
            )

            # Get solutions of substitutions
            sol_eq1 = substituted_eq1_lhs.evalf()
            sol_eq2 = substituted_eq2_lhs.evalf()

            # Compare with exact expected values
            expected_sol_eq1 = her_or_prev
            expected_sol_eq2 = sympy.Integer(1)
            test.assertAlmostEqual(expected_sol_eq1, sol_eq1, accuracy_places)
            test.assertAlmostEqual(expected_sol_eq2, sol_eq2, accuracy_places)
