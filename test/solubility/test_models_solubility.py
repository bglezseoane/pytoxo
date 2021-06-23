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

"""Part of PyToxo solubility test suite."""

import os
import random
import unittest

import pytoxo.calculations
import pytoxo.model
import pytoxo.errors


class ModelsSolubilityTestSuite(unittest.TestCase):
    """This test focuses on verifying that PyToxo is capable of solving at
    least the same models as Toxo. For this, go to the Toxo outputs folders
    (`toxo_outputs/calculate_all_tables_with_times_max_p` and
    `toxo_outputs/calculate_all_tables_with_times_max_h`) and get all
    correct penetrance tables (those with all penetrances between 0 and 1).

    Then, parse the variables (model, MAFs and prevalence or heritability) of
    these cases and assert that here with PyToxo the cases are also solved.

    To be capable to run the test in a reasonable time, as default only some
    cases are randomly selected and used, but editing the flag `exhaustive`
    to true is possible to use all cases.

    To achieve Toxo outputs, some modifications of scripts
    `calculate_all_tables_with_times_max_p.m` and
    `calculate_all_tables_with_times_max_h.m` has been fixed to define cases
    with a larger set of MAFs, prevalences and heritabilities.

    Tests for PyToxo penetrance table generation process at integration
    level.
    """

    # ####################### EDIT HERE #######################
    # Set to true to check all possible cases. Consider that it will take a lot of time
    exhaustive = False
    n_cases_to_check = 20  # If `exhaustive = False`, controls how many cases use
    # #########################################################

    def test_models_solubility(self):
        """This test focuses on verifying that PyToxo is capable of solving at
        least the same models as Toxo. For this, go to the Toxo outputs folders
        (`test/toxo_outputs/calculate_all_tables_with_times_max_p` and
        `test/toxo_outputs/calculate_all_tables_with_times_max_h`) and get all
        correct penetrance tables (those with all penetrances between 0 and 1).

        Then, parse the variables (model, MAFs and prevalence or heritability) of
        these cases and assert that here with PyToxo the cases are also solved.

        To achieve Toxo outputs, some modifications of scripts
        `calculate_all_tables_with_times_max_p.m` and
        `calculate_all_tables_with_times_max_h.m` has been fixed to define cases
        with a larger set of MAFs, prevalences and heritabilities.
        """
        # Get all Toxo output files
        toxo_cases = [
            os.path.join("toxo_outputs", "calculate_all_tables_with_times_max_h", p)
            for p in os.listdir(
                os.path.join("toxo_outputs", "calculate_all_tables_with_times_max_h")
            )
            if p.endswith(".csv")
        ] + [
            os.path.join("toxo_outputs", "calculate_all_tables_with_times_max_p", p)
            for p in os.listdir(
                os.path.join("toxo_outputs", "calculate_all_tables_with_times_max_p")
            )
            if p.endswith(".csv")
        ]

        # Filter corrupted tables
        for toxo_case_path in toxo_cases:
            with open(toxo_case_path, "r") as f:
                toxo_output_content = f.readlines()

            toxo_output_penetrances = [l.split(",")[1] for l in toxo_output_content]
            for toxo_output_penetrance in toxo_output_penetrances:
                if (
                    float(toxo_output_penetrance) < 0
                    or float(toxo_output_penetrance) > 1
                ):
                    toxo_cases.remove(toxo_case_path)
                    break

        """There are too many cases to check due to this repository contains
        archived a lot of Toxo errors. The following lines serve to select a
        limited set of them to run this script. The exhaustive allow a full
        execution is is set to true"""
        if not self.exhaustive:
            toxo_cases = random.choices(
                toxo_cases, k=self.n_cases_to_check
            )  # Select `k` cases randomly
            toxo_cases = sorted(toxo_cases)  # Reorder after selection

        # Print how many cases are going to be checked
        print(f"There going to be run {len(toxo_cases)} solubility case checks.")

        fail_flag = False  # True when at least a case fails
        # Run filtered cases parsing configurations
        for toxo_case_path in toxo_cases:
            toxo_case = os.path.basename(toxo_case_path)
            model = "_".join(toxo_case.split("_")[:2])
            maf = float(toxo_case.split("_")[2])
            rest = toxo_case.split("_")[3].replace(".csv", "")
            if "p" in rest:
                max_method = pytoxo.model.Model.find_max_heritability_table
                her_or_prev = float(rest.replace("p", ""))
            else:
                max_method = pytoxo.model.Model.find_max_prevalence_table
                her_or_prev = float(rest.replace("h", ""))

            # Generate the model
            model = pytoxo.model.Model(os.path.join("models", f"{model}.csv"))

            # Generate penetrance table
            try:
                ptable = max_method(model, [maf] * model.order, her_or_prev)

                # Assert PyToxo achieved table is correct
                for penetrance in ptable._penetrance_values:
                    self.assertGreaterEqual(penetrance, 0)
                    self.assertLessEqual(
                        penetrance, 1.00000000000001
                    )  # TODO: Fix this in a more elegant way
            except pytoxo.errors.ResolutionError:
                # This approach is to run all test also when someone fails
                fail_flag = True
                print(
                    f"[ERROR]: PyToxo cannot solve case {os.path.basename(toxo_case_path)} and Toxo yes."
                )
            except AssertionError:
                # This approach is to run all test also when someone fails
                fail_flag = True
                print(
                    f"[ERROR]: PyToxo table for case {os.path.basename(toxo_case_path)} is corrupted."
                )

        # Finally check fail flag
        self.assertFalse(fail_flag)
