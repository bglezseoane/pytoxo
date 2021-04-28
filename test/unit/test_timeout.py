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

"""PyToxo timeout unit test suite."""

import os
import unittest

import sympy

import pytoxo.calculations
import pytoxo.errors
import pytoxo.model


class ModelUnitTestSuite(unittest.TestCase):
    """Tests for the timeout handling used during the solving process."""

    def test_timeout_raising_solve(self):
        """Test the raising of a timeout exception during the solution of a
        big model with a too tiny timeout limit.

        Warnings: These test are a bit speculative. They use checks based on
        the premise that the operations involved cannot be performed in such
        a short time. These times are very very insignificant compared to the
        necessary ones, but they still make the tests speculative and not
        based on more deterministic logical conditions. Who knows if in ten
        years these operations will not be able to be solved in those
        computational times...
        """
        m = pytoxo.model.Model(os.path.join("models", "multiplicative_8.csv"))

        long_sample_equation = (
            "-625000000000000000000000000000000000000000000000000000000"
            "*(395868833065051*(x*(y + 1) - x*(98967208266262800000000*y + "
            "46221064723759*(y + 1)**256 + 5423271594254400*(y + 1)**128 + "
            "278394608505060000*(y + 1)**64 + 8166241849481750000*(y + 1)**32 "
            "+ 149714433907165000000*(y + 1)**16 + 1756649357844070000000*(y "
            "+ 1)**8 + 12882095290856500000000*(y + 1)**4 + "
            "53982113599779800000000*(y + 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000 + 215928454399119*(x*(y + 1)**2 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000 + 515283811634261*(x*(y + 1)**4 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/1000000000000000000000 + 702659743137629*(x*(y + 1)**8 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/10000000000000000000000 + 598857735628661*(x*(y + 1)**16 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000000 + 32664967397927*(x*(y + 1)**32 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000000 + 13919730425253*(x*(y + 1)**64 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/1250000000000000000000000 + 3389544746409*(x*(y + 1)**128 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/15625000000000000000000000 + 46221064723759*(x*(y + 1)**256 - "
            "x*(98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/25000000000000000000000000000 + 499996645075379*(-x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000 + "
            "x)**2/500000000000000)/(x*(x*(98967208266262800000000*y + "
            "46221064723759*(y + 1)**256 + 5423271594254400*(y + 1)**128 + "
            "278394608505060000*(y + 1)**64 + 8166241849481750000*(y + 1)**32 "
            "+ 149714433907165000000*(y + 1)**16 + 1756649357844070000000*(y "
            "+ 1)**8 + 12882095290856500000000*(y + 1)**4 + "
            "53982113599779800000000*(y + 1)**2 + "
            "24999931220977200000000000000) - "
            "25000000000000000000000000000)*(98967208266262800000000*y + "
            "46221064723759*(y + 1)**256 + 5423271594254400*(y + 1)**128 + "
            "278394608505060000*(y + 1)**64 + 8166241849481750000*(y + 1)**32 "
            "+ 149714433907165000000*(y + 1)**16 + 1756649357844070000000*(y "
            "+ 1)**8 + 12882095290856500000000*(y + 1)**4 + "
            "53982113599779800000000*(y + 1)**2 + "
            "24999931220977200000000000000)) "
        )

        with self.assertRaises(pytoxo.errors.ResolutionError) as e:
            # This process takes more than 30 minutes in a powerful six-core
            # personal computer) tiny timeout
            tiny_timeout = 1

            # m.find_max_prevalence_table([0.12] * 8, 0.95, solve_timeout=tiny_timeout)
            # The above commented line generates an equation similar to
            # `long_sample_equation`, which is hardcoded to avoid this test
            # takes a long time. Uncomment above line and comment following
            # one to an exhaustive execution
            m._solve([sympy.sympify(long_sample_equation)], solve_timeout=tiny_timeout)
        self.assertEqual(e.exception.cause, "Exceeded timeout")

    def test_timeout_aborting_simplification(self):
        """Test that a simplification of the `calculations` module is aborted
        due to the timeout."""
        long_sample_equation_str = (
            "-625000000000000000000000000000000000000000000000000000000"
            "*(395868833065051*(x*(y + 1) - x*(98967208266262800000000*y + "
            "46221064723759*(y + 1)**256 + 5423271594254400*(y + 1)**128 + "
            "278394608505060000*(y + 1)**64 + 8166241849481750000*(y + 1)**32 "
            "+ 149714433907165000000*(y + 1)**16 + 1756649357844070000000*(y "
            "+ 1)**8 + 12882095290856500000000*(y + 1)**4 + "
            "53982113599779800000000*(y + 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000 + 215928454399119*(x*(y + 1)**2 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000 + 515283811634261*(x*(y + 1)**4 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/1000000000000000000000 + 702659743137629*(x*(y + 1)**8 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/10000000000000000000000 + 598857735628661*(x*(y + 1)**16 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000000 + 32664967397927*(x*(y + 1)**32 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/100000000000000000000000 + 13919730425253*(x*(y + 1)**64 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/1250000000000000000000000 + 3389544746409*(x*(y + 1)**128 - x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/15625000000000000000000000 + 46221064723759*(x*(y + 1)**256 - "
            "x*(98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000)**2"
            "/25000000000000000000000000000 + 499996645075379*(-x*("
            "98967208266262800000000*y + 46221064723759*(y + 1)**256 + "
            "5423271594254400*(y + 1)**128 + 278394608505060000*(y + 1)**64 + "
            "8166241849481750000*(y + 1)**32 + 149714433907165000000*(y + "
            "1)**16 + 1756649357844070000000*(y + 1)**8 + "
            "12882095290856500000000*(y + 1)**4 + 53982113599779800000000*(y "
            "+ 1)**2 + "
            "24999931220977200000000000000)/25000000000000000000000000000 + "
            "x)**2/500000000000000)/(x*(x*(98967208266262800000000*y + "
            "46221064723759*(y + 1)**256 + 5423271594254400*(y + 1)**128 + "
            "278394608505060000*(y + 1)**64 + 8166241849481750000*(y + 1)**32 "
            "+ 149714433907165000000*(y + 1)**16 + 1756649357844070000000*(y "
            "+ 1)**8 + 12882095290856500000000*(y + 1)**4 + "
            "53982113599779800000000*(y + 1)**2 + "
            "24999931220977200000000000000) - "
            "25000000000000000000000000000)*(98967208266262800000000*y + "
            "46221064723759*(y + 1)**256 + 5423271594254400*(y + 1)**128 + "
            "278394608505060000*(y + 1)**64 + 8166241849481750000*(y + 1)**32 "
            "+ 149714433907165000000*(y + 1)**16 + 1756649357844070000000*(y "
            "+ 1)**8 + 12882095290856500000000*(y + 1)**4 + "
            "53982113599779800000000*(y + 1)**2 + "
            "24999931220977200000000000000)) "
        )
        # It is known that Sympy is capable to simplify this equation
        long_sample_equation = sympy.sympify(long_sample_equation_str)

        # Now try to simplify with a known as insufficient (this process
        # takes more than 15 minutes in a powerful six-core personal computer)
        # tiny timeout
        tiny_timeout = 1
        returned_equation = pytoxo.calculations._try_to_simplify(
            long_sample_equation, timeout=tiny_timeout
        )

        # Check that the simplification has not been completed
        self.assertEqual(long_sample_equation, returned_equation)
