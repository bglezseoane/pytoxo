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

This test maximizes prevalence.
"""

import datetime
import os
import platform
import sys

import git
import psutil
import tabulate

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


# Read outputs folder
path = os.path.join("test", "toxo_outputs", "calculate_all_tables_with_times_max_p")
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
    heritability = float(
        model_filename.split("_")[3].replace("h", "").replace(".csv", "")
    )

    # Calculate computation time average
    computation_times = []
    with open(os.path.join(path, times_file), "r") as f:
        times_file_content = f.read().splitlines()
    # Filter relative content
    try:
        times_file_content = [l for l in times_file_content if l.startswith(model_name)]
        times_file_content = [l for l in times_file_content if f"_{maf[0]}_" in l]
        times_file_content = [l for l in times_file_content if f"_h{heritability}" in l]
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
    `find_max_heritability_table`. This time, the input is
    formed by numerical values, so the final solution will be
    numerical and not symbolic."""
    recalculated_heritability = pytoxo.calculations.compute_heritability(
        penetrances, maf, model_order=model.order
    )

    # Compare recalculated and initial heritabilitys
    delta = heritability - recalculated_heritability

    # Append to delta list for the automatic checks
    deltas.append(delta)

    # Append results to the table
    table_content.append(
        [
            model_filename.split("_")[0].capitalize(),
            model_order,
            maf[0],
            heritability,
            delta,
            f"{round(computation_time_av, 4)}",
        ]
    )

# Save the generated report
table_headers = [
    "Model",
    "Order",
    "MAF",
    "Heritability",
    "Error",
    f"Time (s) avg.",
]
now = datetime.datetime.now()
# Calculate file name based in current test name and datetime
script_name = str(os.path.basename(__file__)).split(".")[0]
now = f"{now.year}-{now.month}-{now.day}_{now.hour}:{now.minute}:{now.second}"
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
machine_info_tex = machine_info.replace("_", "\\_")
with open(filename, "x") as f:
    """Paste the table inside a basic document template to can easily
    print it as PDF"""
    f.write(
        "\\documentclass{article}\n"
        "\\usepackage{float}\n"
        "\\begin{document}\n"
        "\\section*{PyToxo Test Suite Report}\n"
        f"\\subsection*{{\\texttt{{{script_name}}}}}\n"
        f"Generated report:\n"
        "\\begin{figure}[H]\n"
        "\\centering\n"
        "\n"
        f"{final_table}"
        "\n"
        "\\caption{Accuracies of the the calculated values for the "
        "\\penetrances by Toxo. Maximizing prevalence}\n"
        "\\end{figure}\n"
        f"Datetime: {now_tex}\n\n"
        f"Machine: \\texttt{{{machine_info_tex}}}\n\n"
        f"Git commit hash: \\texttt{{{git_hash}}}\n\n"
        "\\end{document}"
    )
