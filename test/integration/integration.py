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


class PenetranceTableGenerationTestSuite(unittest.TestCase):
    """Tests for PyToxo penetrance table generation process at integration
    level.
    """

    def test_toxo_calculate_tables(self):
        """Test which simulates the same penetrance tables generation in
        PyToxo than the `example/calculate_tables.m` script in Toxo,
        and compare the outputs of the two programs to be equivalent.
        """
        models_to_test = [
            "additive_3.csv",
            "multiplicative_3.csv",
            "threshold_3.csv",
            "additive_4.csv",
            "multiplicative_4.csv",
            "threshold_4.csv",
        ]  # Models folder. Same models that used in Toxo's `example/calculate_tables.m`

        # Some values to run here and with Toxo
        input_mafs = [0.1, 0.4]
        input_heritabilities = [0.1, 0.8]

        """Output collection from Toxo's `example/calculate_tables.m`, tagged
        in the titles"""
        expected_output_files = [
            "test/integration/toxo_outputs/additive_3_0.1_h0.1.txt",
            "test/integration/toxo_outputs/additive_3_0.1_h0.8.txt",
            "test/integration/toxo_outputs/additive_3_0.4_h0.1.txt",
            "test/integration/toxo_outputs/additive_3_0.4_h0.8.txt",
            "test/integration/toxo_outputs/additive_4_0.1_h0.1.txt",
            "test/integration/toxo_outputs/additive_4_0.1_h0.8.txt",
            "test/integration/toxo_outputs/additive_4_0.4_h0.1.txt",
            "test/integration/toxo_outputs/additive_4_0.4_h0.8.txt",
            "test/integration/toxo_outputs/multiplicative_3_0.1_h0.1.txt",
            "test/integration/toxo_outputs/multiplicative_3_0.1_h0.8.txt",
            "test/integration/toxo_outputs/multiplicative_3_0.4_h0.1.txt",
            "test/integration/toxo_outputs/multiplicative_3_0.4_h0.8.txt",
            "test/integration/toxo_outputs/multiplicative_4_0.1_h0.1.txt",
            "test/integration/toxo_outputs/threshold_3_0.1_h0.1.txt",
            "test/integration/toxo_outputs/threshold_3_0.1_h0.8.txt",
            "test/integration/toxo_outputs/threshold_3_0.4_h0.1.txt",
            "test/integration/toxo_outputs/threshold_3_0.4_h0.8.txt",
        ]

        # Generate the PyToxo output for the example cases
        output_files = []
        """Create a temporal directory where save the output during the test 
        execution"""
        with tempfile.TemporaryDirectory() as output_root:
            for model_name in models_to_test:
                # Generate the model
                model = pytoxo.model.Model(f"models/{model_name}")

                for maf in input_mafs:
                    for heritability in input_heritabilities:
                        # Generate penetrance table
                        ptable = model.find_max_prevalence([maf] * 3, heritability)

                        # Save table to a file to compare then with the Toxo equivalent
                        output_file = os.path.join(
                            output_root,
                            f"{model_name.split('.')[0]}_{maf}_h{heritability}",
                        )
                        ptable.write_to_file(output_file)
                        output_files.append(f"{output_file}.csv")

            for expected_output_file, output_file in zip(
                expected_output_files, expected_output_files
            ):
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
                self.assertEqual(expected_genotypes, genotypes)
                # Check the penetrance with a calculation deviation margin
                for expected_penetrance, penetrance in zip(
                    expected_penetrances, penetrances
                ):
                    self.assertAlmostEqual(
                        float(expected_penetrance), float(penetrance)
                    )


if __name__ == "__main__":
    unittest.main()
