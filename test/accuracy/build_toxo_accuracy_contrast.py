#!/usr/bin/env python

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

Script to build a contrast table with the imported Toxo outputs collection.

This test is generalized to easily edit some values and achieve different
reports. It is possible to maximize both prevalence and heritability. Revise
`EDIT HERE` section of the code.
"""

import datetime
import os
import platform
import sys

import git
import psutil
import tabulate

# Workaround to run script inside a project using project
sys.path.insert(0, "../..")
try:
    import pytoxo
    import pytoxo.calculations
    import pytoxo.model
except ImportError as e:
    raise e
# Assert project home directory because sources have path since here
if os.path.basename(os.getcwd()) == "test":
    os.chdir(os.pardir)
elif os.path.basename(os.getcwd()) == "accuracy":
    os.chdir(os.path.join(os.pardir, os.pardir))

# ####################### EDIT HERE #######################
# Comment or uncomment firsts or seconds of each pair
prev_or_her_str = "Prevalence"
# prev_or_her_str = "Heritability"
recalc_method = pytoxo.calculations.compute_heritability
# recalc_method = pytoxo.calculations.compute_prevalence
prev_or_her_letter = "p"
# prev_or_her_letter = "h"
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


# Read outputs folder
path = os.path.join(
    "test", "toxo_outputs", f"calculate_all_tables_with_times_max_{prev_or_her_letter}"
)
files = os.listdir(path)
times_file = "times.txt"
models = sorted([f for f in files if f.endswith(".csv")])

# Latex table report content
table_content = []

deltas = []  # Delta list for automatic checks

for model_filename in models:
    model_name = "_".join(model_filename.split("_")[:2])
    model_order = int(model_filename.split("_")[1])
    maf = [float(model_filename.split("_")[2])] * model_order
    prev_or_her = float(
        model_filename.split("_")[3]
        .replace(prev_or_her_letter_op, "")
        .replace(".csv", "")
    )

    # Calculate computation time average
    computation_times = []
    with open(os.path.join(path, times_file), "r") as f:
        times_file_content = f.read().splitlines()
    # Filter relative content
    try:
        times_file_content = [l for l in times_file_content if l.startswith(model_name)]
        times_file_content = [l for l in times_file_content if f"_{maf[0]}_" in l]
        times_file_content = [
            l
            for l in times_file_content
            if f"_{prev_or_her_letter_op}{prev_or_her}" in l
        ]
    except IndexError:
        times_file_content = []
    # Filter time values
    if times_file_content:
        computation_times.extend([l.split(":")[-1] for l in times_file_content])

    computation_time_av = sum([float(t) for t in computation_times]) / len(
        computation_times
    )

    # Get penetrance values
    with open(os.path.join(path, model_filename), "r") as f:
        model_file_content = f.read().splitlines()
    penetrances = [float(l.split(",")[1]) for l in model_file_content]

    # Generate PyToxo model
    model = pytoxo.model.Model(os.path.join("models", f"{model_name}.csv"))

    """Use the Toxo penetrance values to run a second
    calculation of the left hand side of the first equation
    of the equation system generated in PyToxo's
    `find_max_heritability_table`. or `find_max_prevalence_table`.
    This time, the input is formed by numerical values, 
    so the final solution will be numerical and not symbolic."""
    recalculated_prev_or_her = recalc_method(penetrances, maf, model_order=model.order)

    # Compare recalculated and initial heritabilitys
    delta = abs(prev_or_her - recalculated_prev_or_her)

    # Append to delta list for the automatic checks
    deltas.append(delta)

    # Append results to the table
    table_content.append(
        [
            model_filename.split("_")[0].capitalize(),
            model_order,
            maf[0],
            prev_or_her,
            delta,
            f"{round(computation_time_av, 4)}",
        ]
    )

# Save the generated report
table_headers = [
    "Model",
    "Order",
    "MAF",
    f"{prev_or_her_str_op}",
    "Error",
    f"Time (s) avg.",
]
now = datetime.datetime.now()
# Calculate file name based in current test name and datetime
script_name = str(os.path.basename(__file__)).split(".")[0]
now = (
    f"{now.year:04}-{now.month:02}-{now.day:02}_{now.hour:02}"
    f"-{now.minute:02}-{now.second:02}"
)
filename = os.path.join(
    "test",
    "accuracy",
    "reports",
    f"{script_name}_{now}.tex",
)
final_table = tabulate.tabulate(
    [table_headers] + table_content, headers="firstrow", tablefmt="latex"
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
script_name_tex = script_name.replace("_", "\\_")
now_tex = now.replace("_", "\\_")
final_table_tex = (
    final_table.replace("tabular", "longtable")
    .replace("\\begin{longtable}", "\\begin{longtable}[H]")
    .replace("\\end{longtable}", "")
)
machine_info_tex = machine_info.replace("_", "\\_")
with open(filename, "x") as f:
    """Paste the table inside a basic document template to can easily
    print it as PDF"""
    f.write(
        "\\documentclass{article}\n"
        "\\usepackage{float}\n"
        "\\usepackage{longtable}\n"
        "\\begin{document}\n"
        "\\section*{PyToxo Test Suite Report}\n"
        f"\\subsection*{{\\texttt{{{script_name_tex}}}}}\n"
        f"Generated report:\n"
        "\n"
        f"{final_table_tex}"
        "\n"
        "\\caption{Accuracies of the the calculated values for the "
        f"penetrances by Toxo. Maximizing {prev_or_her_str.lower()}}}\n"
        "\\end{longtable}\n"
        f"Datetime: {now_tex}\n\n"
        f"Machine: \\texttt{{{machine_info_tex}}}\n\n"
        f"Git commit hash: \\texttt{{{git_hash}}}\n\n"
        "\\end{document}"
    )
