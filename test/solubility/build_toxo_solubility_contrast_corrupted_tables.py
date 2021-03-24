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

"""Part of PyToxo solubility test suite.

Script to build a contrast table with the imported Toxo outputs collection.

Get the Toxo corrupted tables and analyze which cases are solvable by
PyToxo, composing a table.
"""

import datetime
import os
import platform
import sys

import git
import psutil
import tabulate

# Assert project home directory because sources have path since here
if os.path.basename(os.getcwd()) == "test":
    os.chdir(os.pardir)
elif os.path.basename(os.getcwd()) == "solubility":
    os.chdir(os.path.join(os.pardir, os.pardir))

# Workaround to run script inside a project using project
try:
    import pytoxo
    import pytoxo.calculations
    import pytoxo.errors
    import pytoxo.model
except ImportError or ModuleNotFoundError:
    try:
        sys.path.insert(0, os.getcwd())  # Already updated to project root in above step
        import pytoxo
        import pytoxo.calculations
        import pytoxo.errors
        import pytoxo.model
    except ImportError or ModuleNotFoundError as e:
        raise e

# ####################### EDIT HERE #######################
# Comment or uncomment firsts or seconds of each pair
prev_or_her_str = "Prevalence"
# prev_or_her_str = "Heritability"
calc_method = pytoxo.model.Model.find_max_prevalence_table
# calc_method = pytoxo.model.Model.find_max_heritability_table
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


# Some definitions about the Toxo outputs environment
path = os.path.join(
    "test", "toxo_outputs", f"calculate_all_tables_with_times_max_{prev_or_her_letter}"
)
files = os.listdir(path)
outputs = sorted([f for f in files if f.endswith(".csv")])

# Filter models: get only models bad computed by Toxo (corrupted penetrances)
outputs_to_check = []
for output in outputs:
    with open(os.path.join(path, output), "r") as f:
        output_content = f.readlines()
    output_penetrances = [l.split(",")[1] for l in output_content]
    for output_penetrance in output_penetrances:
        if float(output_penetrance) > 1 or float(output_penetrance) < 0:
            outputs_to_check.append(output)

# Latex table report content
table_content = []

deltas = []  # Delta list for automatic checks

try:
    for output in outputs_to_check:
        model_name = "_".join(output.split("_")[:2])
        model_order = int(output.split("_")[1])
        maf = [float(output.split("_")[2])] * model_order
        prev_or_her = float(
            output.split("_")[3].replace(prev_or_her_letter_op, "").replace(".csv", "")
        )

        # Generate PyToxo model
        model = pytoxo.model.Model(os.path.join("models", f"{model_name}.csv"))

        # Try to solve with PyToxo and annotate
        try:
            pt = calc_method(model, maf, prev_or_her, check=True)

            # Append results to the table if PyToxo finish on success the attempt
            table_content.append(
                [
                    model_name.capitalize(),
                    model_order,
                    maf[0],
                    prev_or_her,
                ]
            )
        except pytoxo.errors.ResolutionError or pytoxo.errors.UnsolvableModelError:
            continue  # Next case...
except KeyboardInterrupt:
    pass  # Continue to save a report with the already calculated data

# Check if there at least a case to list
if not table_content:
    print("There are not any case insolvable with Toxo and yes with PyToxo.")
    sys.exit(0)

# Save the generated report
table_headers = [
    "Model",
    "Order",
    "MAF",
    f"{prev_or_her_str_op}",
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
    "solubility",
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
        "\\caption{List of models that Toxo cannot solve (calculates a "
        f"corrupted table) but PyToxo yes. Maximizing {prev_or_her_str.lower()}}}\n"
        "\\end{longtable}\n"
        f"Datetime: {now_tex}\n\n"
        f"Machine: \\texttt{{{machine_info_tex}}}\n\n"
        f"Git commit hash: \\texttt{{{git_hash}}}\n\n"
        "\\end{document}"
    )
