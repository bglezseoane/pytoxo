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

"""PyToxo integration test suite."""

import os
import tempfile
import unittest

import pytoxo.calculations
import pytoxo.model


class ToxoCalculateTablesTestSuite(unittest.TestCase):
    """Test suite which simulates the same penetrance tables generation in
    PyToxo than the `example/calculate_tables.m` script in Toxo, and compare
    the outputs of the two programs to be equivalent.

    Tests for PyToxo penetrance table generation process at integration
    level.
    """

    def test_toxo_calculate_tables_additive_3_maf_1_h_1(self):
        """Test the calculation of a penetrance table for the `additive_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_3"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_3_maf_1_h_1(self):
        """Test the calculation of a penetrance table for the `multiplicative_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """

        model = "multiplicative_3"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_3_maf_1_h_1(self):
        """Test the calculation of a penetrance table for the `threshold_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_3"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_4_maf_1_h_1(self):
        """Test the calculation of a penetrance table for the `additive_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_4"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_4_maf_1_h_1(self):
        """Test the calculation of a penetrance table for the `multiplicative_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "multiplicative_4"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_4_maf_1_h_1(self):
        """Test the calculation of a penetrance table for the `threshold_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_4"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_3_maf_1_h_8(self):
        """Test the calculation of a penetrance table for the `additive_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_3"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_3_maf_1_h_8(self):
        """Test the calculation of a penetrance table for the `multiplicative_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """

        model = "multiplicative_3"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_3_maf_1_h_8(self):
        """Test the calculation of a penetrance table for the `threshold_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_3"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_4_maf_1_h_8(self):
        """Test the calculation of a penetrance table for the `additive_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_4"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_4_maf_1_h_8(self):
        """Test the calculation of a penetrance table for the `multiplicative_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "multiplicative_4"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_4_maf_1_h_8(self):
        """Test the calculation of a penetrance table for the `threshold_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_4"
        order = int(model.split("_")[1])
        maf = 0.1
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_3_maf_4_h_1(self):
        """Test the calculation of a penetrance table for the `additive_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_3"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_3_maf_4_h_1(self):
        """Test the calculation of a penetrance table for the `multiplicative_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """

        model = "multiplicative_3"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_3_maf_4_h_1(self):
        """Test the calculation of a penetrance table for the `threshold_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_3"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_4_maf_4_h_1(self):
        """Test the calculation of a penetrance table for the `additive_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_4"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_4_maf_4_h_1(self):
        """Test the calculation of a penetrance table for the `multiplicative_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "multiplicative_4"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_4_maf_4_h_1(self):
        """Test the calculation of a penetrance table for the `threshold_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_4"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.1
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_3_maf_4_h_8(self):
        """Test the calculation of a penetrance table for the `additive_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_3"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_3_maf_4_h_8(self):
        """Test the calculation of a penetrance table for the `multiplicative_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """

        model = "multiplicative_3"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_3_maf_4_h_8(self):
        """Test the calculation of a penetrance table for the `threshold_3`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_3"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_additive_4_maf_4_h_8(self):
        """Test the calculation of a penetrance table for the `additive_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "additive_4"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_multiplicative_4_maf_4_h_8(self):
        """Test the calculation of a penetrance table for the `multiplicative_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "multiplicative_4"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    def test_toxo_calculate_tables_threshold_4_maf_4_h_8(self):
        """Test the calculation of a penetrance table for the `threshold_4`
        model, with a MAF 0.1 and an heritability 0.1. To check the result,
        compares it with the collected Toxo's output for the same input.
        """
        model = "threshold_4"
        order = int(model.split("_")[1])
        maf = 0.4
        heritability = 0.8
        expected_output_file = (
            f"test/integration/toxo_outputs/{model}_{maf}_h{heritability}.txt"
        )

        # Run within the helper function
        self.helper_toxo_calculate_tables(
            test=self,
            model_file=f"models/{model}.csv",
            expected_output_file=expected_output_file,
            maf=maf,
            heritability=heritability,
            accuracy_delta=0.01,
        )

    @staticmethod
    def helper_toxo_calculate_tables(
        test, model_file, expected_output_file, maf, heritability, accuracy_delta
    ):
        """Helper method with the test skeleton to the penetrance tables
        generation tests of `ToxoCalculateTablesTestSuite` test suite.
        """

        """Create a temporal directory where save the output during the test 
        execution"""
        with tempfile.TemporaryDirectory() as output_root:
            # Generate the model
            model = pytoxo.model.Model(model_file)

            # Generate penetrance table
            ptable = model.find_max_prevalence([maf] * model.order, heritability)

            # Save table to a file to compare then with the Toxo equivalent
            output_file = os.path.join(
                output_root,
                f"{os.path.basename(model_file)}_{maf}_h{heritability}.txt",
            )
            ptable.write_to_file(output_file)

            # Compare Toxo and PyToxo outputs
            with open(expected_output_file, "r") as f:
                expected_output = f.readlines()

            # Split genotypes and penetrance values columns
            expected_genotypes = [line.split(",")[0] for line in expected_output]
            expected_penetrances = [line.split(",")[1] for line in expected_output]

            with open(output_file, "r") as f:
                output = f.readlines()

            # Split genotypes and penetrance values columns
            genotypes = [line.split(",")[0] for line in output]
            penetrances = [line.split(",")[1] for line in output]

            # Genotypes should be exactly the same
            test.assertEqual(expected_genotypes, genotypes)
            # Check the penetrance with a calculation deviation margin
            for expected_penetrance, penetrance in zip(
                expected_penetrances, penetrances
            ):
                # Cast the penetrances
                expected_penetrance = float(expected_penetrance)
                penetrance = float(penetrance)
                # First check penetrance values are coherent
                test.assertGreaterEqual(1, expected_penetrance)
                test.assertLessEqual(0, expected_penetrance)
                test.assertGreaterEqual(1, penetrance)
                test.assertLessEqual(0, penetrance)
                # Calculate delta between the two penetrances
                penetrance_delta = abs(expected_penetrance - penetrance)
                test.assertLess(penetrance_delta, accuracy_delta)


if __name__ == "__main__":
    unittest.main()
