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

"""Part of the PyToxo accuracy test suite.

Script to check PyToxo calculated penetrances accuracy and computation
time for different models.

These scripts run some automatic checks and generate tables as Latex format to
revise manually the accuracies and times over different PyToxo versions.
"""

import datetime
import os
import platform
import time
import unittest

import psutil

import pytoxo
import pytoxo.calculations
import pytoxo.model
import tabulate

_TOLERABLE_ACCURACY_DELTA = 1e-15  # Maximum error
_TEST_REPETITIONS = 1  # To confirm computation times with an average


class PenetrancesAccuracyMaxPTestSuite(unittest.TestCase):
    """Test suite to check PyToxo calculated penetrances accuracy and
    computation time for different models. I.e. the final penetrances of the
    penetrances tables.

    This test maximizes prevalence.
    """

    def test_penetrances_accuracy(self):
        models = [
            # "additive_2",
            "additive_3",
            "additive_4",
            # "additive_5",
            # "additive_6",
            # "additive_7",
            # "additive_8",
            # "multiplicative_2",
            "multiplicative_3",
            "multiplicative_4",
            # "multiplicative_5",
            # "multiplicative_6",
            # "multiplicative_7",
            # "multiplicative_8",
            # "threshold_2",
            "threshold_3",
            "threshold_4",
            # "threshold_5",
            # "threshold_6",
            # "threshold_7",
            # "threshold_8",
        ]  # Uncomment ones to use in the test
        mafs = [0.1, 0.4]
        heritabilities = [0.1, 0.8]

        # Latex table report content
        table_headers = [
            "Model",
            "Order",
            "MAF",
            "Heritability",
            "Error",
            f"Time (s) avg. {_TEST_REPETITIONS}",
        ]
        table_content = [table_headers]

        deltas = []  # Delta list for automatic checks

        for model_name in models:
            model_order = int(model_name.split("_")[1])
            for maf in mafs:
                maf = [maf] * model_order
                for heritability in heritabilities:
                    # Generate model
                    model = pytoxo.model.Model(
                        os.path.join("models", f"{model_name}.csv")
                    )

                    # Run the computation configured times, saving consumed times
                    computation_times = []
                    for _ in range(_TEST_REPETITIONS):
                        t0 = time.time()
                        # Full penetrance table process since model generation
                        ptable_sol = model.find_max_prevalence_table(maf, heritability)
                        tf = time.time()
                        computation_times.append(tf - t0)

                    # Calculate computation time average
                    computation_time_av = sum(computation_times) / len(
                        computation_times
                    )

                    """Use the obtained penetrance values to run a second 
                    calculation of the left hand side of the first equation 
                    of the equation system generated in
                    `find_max_prevalence_table`. This time, the input is 
                    formed by numerical values, so the final solution will be
                    numerical and not symbolic."""
                    recalculated_heritability = (
                        pytoxo.calculations.compute_heritability(
                            ptable_sol._penetrance_values, maf, model_order=model.order
                        )
                    )

                    # Compare recalculated and initial prevalences
                    delta = heritability - recalculated_heritability

                    # Append to delta list for the automatic checks
                    deltas.append(delta)

                    # Append results to the table
                    table_content.append(
                        [
                            model_name.split("_")[0].capitalize(),
                            model_order,
                            maf[0],
                            heritability,
                            delta,
                            f"{round(computation_time_av, 4)}",
                        ]
                    )

        # Save the generated report
        now = datetime.datetime.now()
        # Calculate file name based in current test name and datetime
        test_name = str(self).split(" ")[0]
        module_name = str(self.__module__)
        now = f"{now.year}-{now.month}-{now.day}_{now.hour}:{now.minute}:{now.second}"
        filename = os.path.join(
            "test",
            "accuracy",
            "reports",
            f"{test_name}_{now}.tex",
        )
        final_table = tabulate.tabulate(
            table_content, headers="firstrow", tablefmt="latex"
        )
        machine_info = (
            f"{platform.platform()}, "
            f"{psutil.cpu_count(logical=True)} core, "
            f"{psutil.cpu_count(logical=False)} physical core, "
            f"{psutil.cpu_freq().max:.2f} MHz max freq."
        )
        # Some Latex patches
        test_name_tex = test_name.replace("_", "\\_")
        module_name_tex = module_name.replace("_", "\\_")
        now_tex = now.replace("_", "\\_")
        machine_info_tex = machine_info.replace("_", "\\_")
        with open(filename, "x") as f:
            """Paste the table inside a basic document template to can easily
            print it as PDF"""
            f.write(
                "\\documentclass{article}\n"
                "\\usepackage{float}\n"
                "\\begin{document}\n"
                "\\section*{PyToxo Test Suite Report}\n"
                f"\\subsection*{{\\texttt{{{module_name_tex}: {test_name_tex}}}}}\n"
                f"Generated report:\n"
                "\\begin{figure}[H]\n"
                "\\centering\n"
                "\n"
                f"{final_table}"
                "\n"
                "\\caption{Accuracies of the the calculated values for the "
                "\\texttt{x} and \\texttt{y} variables used in the penetrance "
                "tables}\n"
                "\\end{figure}\n"
                f"Datetime: {now_tex}\n\n"
                f"Machine: \\texttt{{{machine_info_tex}}}\n"
                "\\end{document}"
            )

        # Automatic checks against configured tolerable delta
        for delta in deltas:
            self.assertGreaterEqual(_TOLERABLE_ACCURACY_DELTA, delta)
