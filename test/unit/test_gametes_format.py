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
import unittest

import pytoxo.errors
import pytoxo.model


class GAMETESFormatTestSuite(unittest.TestCase):
    """Tests for check to correct composition of the `PTable` formatted as
    GAMETES format.
    """

    def test_ptable_as_gametes_check_disposition_as_unknown_2(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a valid sample file to the same input. Only compares
        the table disposition, not the values due to the sample table is not
        generated with a known repository model file."""
        test_order = 2  # Determine the sources to use in this particular test

        sample_filename = os.path.join(
            "test", "unit", "gametes_output_samples", f"unknown_{test_order}.txt"
        )
        with open(sample_filename, "r") as sample_file:
            sample = sample_file.readlines()

        # Parse the file to compose the experiment
        for l in sample:
            if l.startswith("Minor allele frequencies:\t"):
                mafs = [
                    float(maf)
                    for maf in l.replace("Minor allele frequencies:\t", "").split("\t")
                ]
                break
        for l in sample:
            if l.startswith("Heritability: "):
                h = float(l.replace("Heritability: ", ""))
                break
        for l in sample:
            if l.startswith("Table:"):
                expected_output_table = sample[sample.index(l) + 2 :]
                break

        # Compound the model with the read confifuration and generate the table
        m = pytoxo.model.Model(
            os.path.join(
                "models", f"additive_{test_order}.csv"
            )  # Unknown so unchecked, only use any one
        )
        pt = m.find_max_prevalence_table(mafs=mafs, h=h, check=False)
        pt_as_gametes = pt._compound_table_as_gametes().splitlines(keepends=True)

        # Discard headers to compare only the table
        for l in pt_as_gametes:
            if l.startswith("Table:"):
                output_table = pt_as_gametes[pt_as_gametes.index(l) + 2 :]
                break

        # Check only the disposition of the table members
        for expected_output_table_line, output_table_line in zip(
            expected_output_table, output_table
        ):
            if expected_output_table_line == "\n":
                self.assertEqual(expected_output_table_line, output_table_line)
            else:
                for ev, v in zip(
                    expected_output_table_line.split(","), output_table_line.split(",")
                ):
                    self.assertEqual(type(ev), type(v))

    def test_ptable_as_gametes_check_disposition_as_unknown_3(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a valid sample file to the same input. Only compares
        the table disposition, not the values due to the sample table is not
        generated with a known repository model file."""
        test_order = 3  # Determine the sources to use in this particular test

        sample_filename = os.path.join(
            "test", "unit", "gametes_output_samples", f"unknown_{test_order}.txt"
        )
        with open(sample_filename, "r") as sample_file:
            sample = sample_file.readlines()

        # Parse the file to compose the experiment
        for l in sample:
            if l.startswith("Minor allele frequencies:\t"):
                mafs = [
                    float(maf)
                    for maf in l.replace("Minor allele frequencies:\t", "").split("\t")
                ]
                break
        for l in sample:
            if l.startswith("Heritability: "):
                h = float(l.replace("Heritability: ", ""))
                break
        for l in sample:
            if l.startswith("Table:"):
                expected_output_table = sample[sample.index(l) + 2 :]
                break

        # Compound the model with the read confifuration and generate the table
        m = pytoxo.model.Model(
            os.path.join(
                "models", f"additive_{test_order}.csv"
            )  # Unknown so unchecked, only use any one
        )
        pt = m.find_max_prevalence_table(mafs=mafs, h=h, check=False)
        pt_as_gametes = pt._compound_table_as_gametes().splitlines(keepends=True)

        # Discard headers to compare only the table
        for l in pt_as_gametes:
            if l.startswith("Table:"):
                output_table = pt_as_gametes[pt_as_gametes.index(l) + 2 :]
                break

        # Check only the disposition of the table members
        for expected_output_table_line, output_table_line in zip(
            expected_output_table, output_table
        ):
            if expected_output_table_line == "\n":
                self.assertEqual(expected_output_table_line, output_table_line)
            else:
                for ev, v in zip(
                    expected_output_table_line.split(","), output_table_line.split(",")
                ):
                    self.assertEqual(type(ev), type(v))

    def test_ptable_as_gametes_check_disposition_as_unknown_4(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a valid sample file to the same input. Only compares
        the table disposition, not the values due to the sample table is not
        generated with a known repository model file."""
        test_order = 4  # Determine the sources to use in this particular test

        sample_filename = os.path.join(
            "test", "unit", "gametes_output_samples", f"unknown_{test_order}.txt"
        )
        with open(sample_filename, "r") as sample_file:
            sample = sample_file.readlines()

        # Parse the file to compose the experiment
        for l in sample:
            if l.startswith("Minor allele frequencies:\t"):
                mafs = [
                    float(maf)
                    for maf in l.replace("Minor allele frequencies:\t", "").split("\t")
                ]
                break
        for l in sample:
            if l.startswith("Heritability: "):
                h = float(l.replace("Heritability: ", ""))
                break
        for l in sample:
            if l.startswith("Table:"):
                expected_output_table = sample[sample.index(l) + 2 :]
                break

        # Compound the model with the read confifuration and generate the table
        m = pytoxo.model.Model(
            os.path.join(
                "models", f"additive_{test_order}.csv"
            )  # Unknown so unchecked, only use any one
        )
        pt = m.find_max_prevalence_table(mafs=mafs, h=h, check=False)
        pt_as_gametes = pt._compound_table_as_gametes().splitlines(keepends=True)

        # Discard headers to compare only the table
        for l in pt_as_gametes:
            if l.startswith("Table:"):
                output_table = pt_as_gametes[pt_as_gametes.index(l) + 2 :]
                break

        # Check only the disposition of the table members
        for expected_output_table_line, output_table_line in zip(
            expected_output_table, output_table
        ):
            if expected_output_table_line == "\n":
                self.assertEqual(expected_output_table_line, output_table_line)
            else:
                for ev, v in zip(
                    expected_output_table_line.split(","), output_table_line.split(",")
                ):
                    self.assertEqual(type(ev), type(v))

    def test_ptable_as_gametes_check_all_as_toxo_2(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 2)

    def test_ptable_as_gametes_check_all_as_toxo_2(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 2)

    def test_ptable_as_gametes_check_all_as_toxo_3(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 3)

    def test_ptable_as_gametes_check_all_as_toxo_4(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 4)

    def test_ptable_as_gametes_check_all_as_toxo_5(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 5)

    def test_ptable_as_gametes_check_all_as_toxo_6(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 6)

    def test_ptable_as_gametes_check_all_as_toxo_7(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 7)

    def test_ptable_as_gametes_check_all_as_toxo_8(self):
        """Test composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        self._helper_ptable_as_gametes_check_all_as_toxo(self, 8)

    @staticmethod
    def _helper_ptable_as_gametes_check_all_as_toxo(test, test_order):
        """Helper method with the test skeleton for the test of
        the composition of the `PTable` formatted as GAMETES format,
        comparing with a Toxo generated file to the same input. Compares both
        the table disposition and the values. More exhaustive than
        `test_ptable_as_gametes_check_disposition_as_unknown_*` test."""
        sample_filename = os.path.join(
            "test", "unit", "gametes_output_samples", f"toxo_{test_order}.txt"
        )
        with open(sample_filename, "r") as sample_file:
            sample = sample_file.readlines()

        # Parse the file to compose the experiment ans compare then
        for l in sample:
            if l.startswith("Minor allele frequencies:\t"):
                expected_output_mafs = [
                    float(maf)
                    for maf in l.replace("Minor allele frequencies:\t", "").split("\t")
                ]
                break
        for l in sample:
            if l.startswith("Heritability: "):
                expected_output_h = float(l.replace("Heritability: ", ""))
                break
        for l in sample:
            if l.startswith("Prevalence: "):
                expected_output_p = float(l.replace("Prevalence: ", ""))
                break
        for l in sample:
            if l.startswith("x: "):
                expected_output_x = float(l.replace("x: ", ""))
                break
        for l in sample:
            if l.startswith("y: "):
                expected_output_y = float(l.replace("y: ", ""))
                break
        for l in sample:
            if l.startswith("Table:"):
                expected_output_table = sample[sample.index(l) + 2 :]
                break

        # Compound the model with the read configuration and generate the table
        m = pytoxo.model.Model(
            os.path.join("models", f"additive_{test_order}.csv")  # Known
        )
        pt = m.find_max_heritability_table(
            mafs=expected_output_mafs, p=expected_output_p, check=False
        )
        pt_as_gametes = pt._compound_table_as_gametes().splitlines(keepends=True)

        # Parse the output table
        for l in pt_as_gametes:
            if l.startswith("Minor allele frequencies:\t"):
                output_mafs = [
                    float(maf)
                    for maf in l.replace("Minor allele frequencies:\t", "").split("\t")
                ]
                break
        for l in pt_as_gametes:
            if l.startswith("Heritability: "):
                output_h = float(l.replace("Heritability: ", ""))
                break
        for l in pt_as_gametes:
            if l.startswith("Prevalence: "):
                output_p = float(l.replace("Prevalence: ", ""))
                break
        for l in pt_as_gametes:
            if l.startswith("x: "):
                output_x = float(l.replace("x: ", ""))
                break
        for l in pt_as_gametes:
            if l.startswith("y: "):
                output_y = float(l.replace("y: ", ""))
                break
        for l in pt_as_gametes:
            if l.startswith("Table:"):
                output_table = sample[sample.index(l) + 2 :]
                break

        # Check the headers
        test.assertEqual(expected_output_mafs, output_mafs)
        test.assertAlmostEqual(
            expected_output_h,
            output_h,
            places=5,  # Loose, accuracy is not checked here
        )
        test.assertAlmostEqual(
            expected_output_p,
            output_p,
            places=5,  # Loose, accuracy is not checked here
        )
        test.assertAlmostEqual(
            expected_output_x,
            output_x,
            places=5,  # Loose, accuracy is not checked here
        )
        test.assertAlmostEqual(
            expected_output_y,
            output_y,
            places=5,  # Loose, accuracy is not checked here
        )

        # Check the table disposition and values
        for expected_output_table_line, output_table_line in zip(
            expected_output_table, output_table
        ):
            if expected_output_table_line == "\n":
                test.assertEqual(expected_output_table_line, output_table_line)
            else:
                for ev, v in zip(
                    expected_output_table_line.split(","), output_table_line.split(",")
                ):
                    test.assertAlmostEqual(float(ev), float(v))
