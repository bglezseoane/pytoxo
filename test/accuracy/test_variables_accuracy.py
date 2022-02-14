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

"""Part of the PyToxo accuracy test suite.

Script to check PyToxo calculated variable solutions accuracy and computation
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
import tabulate

import pytoxo
import pytoxo.errors
import pytoxo.model

# ####################### EDIT HERE #######################
# Flag to control when the generated reports should be saved
print_reports = False

# Comment or uncomment firsts or seconds of each pair
prev_or_her_str = "Prevalence"
# prev_or_her_str = "Heritability"
recalc_method = pytoxo.model.Model._build_max_prevalence_system
# recalc_method = pytoxo.model.Model._build_max_heritability_system
table_method = pytoxo.model.Model.find_max_prevalence_table
# table_method = pytoxo.model.Model.find_max_heritability_table
prev_or_her_letter = "p"
# prev_or_her_letter = "h"

MAFS = [0.1, 0.4]
# MAFS = [0.1, 0.2, 0.3, 0.4, 0.5]
PREVS_OR_HERS = [0.1, 0.8]
# PREVS_OR_HERS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# Flag to select between LaTeX or CSV report type
# report_extension = ".tex"
report_extension = ".csv"
# #########################################################

# Calculate some necessary opposites
if prev_or_her_str == "Prevalence":
    prev_or_her_str_op = "Heritability"
else:
    prev_or_her_str_op = "Prevalence"
if prev_or_her_letter == "p":
    prev_or_her_letter_op = "h"
else:
    prev_or_her_letter_op = "p"


if print_reports:
    import git

_TEST_REPETITIONS = 3  # To confirm computation times with an average


class VariablesAccuracyTestSuite(unittest.TestCase):
    """Test suite to check PyToxo calculated variable solutions accuracy and
    computation time for different models. I.e. the calculated values for the
    `x` and `y` variables used in the penetrance tables.

    This test maximizes prevalence.
    """

    def test_variables_accuracy(self):
        models = [
            "additive_2",
            "additive_3",
            "additive_4",
            "additive_5",
            "additive_6",
            "additive_7",
            "additive_8",
            "multiplicative_2",
            "multiplicative_3",
            "multiplicative_4",
            "threshold_2",
            "threshold_3",
            "threshold_4",
            "threshold_5",
            "threshold_6",
            "threshold_7",
            "threshold_8",
        ]  # Uncomment ones to use in the test
        mafs = MAFS
        prevs_or_hers = PREVS_OR_HERS

        # Latex table report content
        table_headers = [
            "Model",
            "Order",
            "MAF",
            f"{prev_or_her_str_op}",
            "Error",
            f"Time (s) avg. {_TEST_REPETITIONS}",
        ]
        table_content = [table_headers]

        deltas = []  # `(tolerable_delta, delta)` list for automatic checks

        for model_name in models:
            model_order = int(model_name.split("_")[1])
            for maf in mafs:
                maf = [maf] * model_order
                for prev_or_her in prevs_or_hers:
                    # Generate model
                    model = pytoxo.model.Model(
                        os.path.join("models", f"{model_name}.csv")
                    )

                    # Generate equation system
                    eq_system = recalc_method(model, maf, prev_or_her)

                    # Get the equation system PyToxo solution
                    try:
                        vars_sol = model._solve(eq_system, solve_timeout=False)
                    except pytoxo.errors.UnsolvableModelError or pytoxo.errors.ResolutionError:
                        """If resolution tentative fails, simple go to next
                        case. This test has not the responsibility to check
                        model solubility, only accuracy."""
                        continue

                    """Build the table a specified number of 
                    repetitions, to time the calculation time. This step is 
                    useful only to measure the time of the complete process of 
                    solving the equations and to construct the final table, with 
                    final substitution solution. To measure the precision use 
                    `vars_sol` above, which is not so suitable for measuring 
                    times because it is only one of the isolated steps of the 
                    full process"""
                    computation_times = []
                    try:
                        for _ in range(_TEST_REPETITIONS):
                            t0 = time.time()
                            _ = table_method(
                                model,
                                maf,
                                prev_or_her,
                                check=False,
                                solve_timeout=False,
                            )
                            tf = time.time()
                            computation_times.append(tf - t0)
                    except pytoxo.errors.UnsolvableModelError or pytoxo.errors.ResolutionError:
                        """If resolution tentative fails, simple go to next
                        case. This test has not the responsibility to check
                        model solubility, only accuracy."""
                        continue

                    # Calculate computation time average
                    computation_time_av = sum(computation_times) / len(
                        computation_times
                    )

                    # Recuperate prev. or her. equation unequaled
                    eq1_lhs = eq_system[0].lhs  # Left hand side

                    # Now substitute the calculated solutions in the equation
                    substituted_eq1_lhs = eq1_lhs.subs(
                        {
                            model.variables[0]: vars_sol[model.variables[0]],
                            model.variables[1]: vars_sol[model.variables[1]],
                        }
                    )

                    # Get solution of substitution
                    sol_eq1 = substituted_eq1_lhs.evalf()

                    # Compare with exact expected values calculating a delta
                    delta = abs(prev_or_her - sol_eq1)

                    """Append to the list the tolerable delta for the current 
                    model and the achieved delta"""
                    deltas.append(
                        (model.calculate_tolerable_solution_error_delta(), delta)
                    )

                    # Append results to the table
                    table_content.append(
                        [
                            model_name.split("_")[0].capitalize(),
                            model_order,
                            maf[0],
                            prev_or_her,
                            delta,
                            f"{round(computation_time_av, 4)}",
                        ]
                    )

        if print_reports:
            # Save the generated report
            now = datetime.datetime.now()
            # Calculate file name based in current test name and datetime
            test_name = str(self).split(" ")[0]
            module_name = str(self.__module__)
            now = (
                f"{now.year:04}-{now.month:02}-{now.day:02}_{now.hour:02}"
                f"-{now.minute:02}-{now.second:02}"
            )
            filename = os.path.join(
                "test",
                "accuracy",
                "reports",
                f"{test_name}_max_{prev_or_her_letter}_{now}{report_extension}",
            )
            if report_extension == ".tex":
                final_table = tabulate.tabulate(
                    table_content, headers="firstrow", tablefmt="latex"
                )
                machine_info = (
                    f"{platform.platform()}, "
                    f"{psutil.cpu_count(logical=True)} core, "
                    f"{psutil.cpu_count(logical=False)} physical core, "
                    f"{psutil.cpu_freq().max:.2f} MHz max freq."
                )
                """Retrieve current repository commit reference to locate the report 
                in the history"""
                git_hash = git.Repo(search_parent_directories=True).head.object.hexsha
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
                        "tables. Corrupted tables are discarded. "
                        f"Maximizing {prev_or_her_str.lower()}}}\n"
                        "\\end{figure}\n"
                        f"Datetime: {now_tex}\n\n"
                        f"Machine: \\texttt{{{machine_info_tex}}}\n\n"
                        f"Git commit hash: \\texttt{{{git_hash}}}\n\n"
                        "\\end{document}"
                    )
            else:
                with open(filename, "x") as f:
                    for line in table_content:
                        f.write(";".join([str(e) for e in line]))
                        f.write("\n")

        # Automatic checks against configured tolerable delta
        for tolerable_delta, delta in deltas:
            self.assertGreaterEqual(tolerable_delta, delta)
