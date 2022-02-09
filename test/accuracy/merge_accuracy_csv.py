#!/usr/bin/env python

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

Script to merge CSV outputs of scripts `test_variables_accuracy.py` and
`build_toxo_accuracy_contrast.py` to achieve contrast data between PyToxo and
Toxo.
"""

import datetime
import os
import sys

import numpy as np
import pandas as pd  # TODO: Add to requirements
import seaborn as sns
import matplotlib.pyplot as plt

# ####################### EDIT HERE #######################
PYTOXO_REPORT_FILENAME = "test_variables_accuracy_2022-02-08_03-31-41.csv"
TOXO_REPORT_FILENAME = "build_toxo_accuracy_contrast_2022-02-08_17-51-47.csv"
SAVE = True
PLOT = True
# #########################################################

# Assert project home directory because sources have path since here
if os.path.basename(os.getcwd()) == "test":
    os.chdir(os.pardir)
elif os.path.basename(os.getcwd()) == "accuracy":
    os.chdir(os.path.join(os.pardir, os.pardir))

# Workaround to run script inside a project using project
try:
    import pytoxo
    import pytoxo.calculations
    import pytoxo.model
except ImportError or ModuleNotFoundError:
    try:
        sys.path.insert(0, os.getcwd())  # Already updated to project root in above step
        import pytoxo
        import pytoxo.calculations
        import pytoxo.model
    except ImportError or ModuleNotFoundError as e:
        raise e

REPORTS_PATH = os.path.join("test", "accuracy", "reports")

# Read reports
pytoxo_report = pd.read_csv(
    os.path.join(REPORTS_PATH, PYTOXO_REPORT_FILENAME), delimiter=";"
)
toxo_report = pd.read_csv(
    os.path.join(REPORTS_PATH, TOXO_REPORT_FILENAME), delimiter=";"
)

# Merge reports
merged = pytoxo_report.merge(toxo_report, on=["Model", "Order", "MAF", "Heritability"])

# Calculate file name based in current test name and datetime
now = datetime.datetime.now()
script_name = str(os.path.basename(__file__)).split(".")[0]
now = (
    f"{now.year:04}-{now.month:02}-{now.day:02}_{now.hour:02}"
    f"-{now.minute:02}-{now.second:02}"
)
output_report_path = os.path.join(
    "test",
    "accuracy",
    "reports",
    f"{script_name}_{now}.csv",
)

# Write final file
if SAVE:
    with open(output_report_path, "w") as f:  # Patch to set custom headers
        f.write(
            ";".join(
                [
                    "Model",
                    "Order",
                    "MAF",
                    "Heritability",
                    "Error PyToxo",
                    "Time (s) avg. 3 PyToxo",
                    "Error Toxo",
                    "Time (s) avg. 3 Toxo\n",
                ]
            )
        )
    merged.to_csv(
        output_report_path,
        header=False,
        index=False,
        sep=";",
        mode="a",
    )

if PLOT:
    # Calculate global accumulated errors
    pytoxo_errors = merged.Error_x
    toxo_errors = merged.Error_y
    # toxo_error = toxo_error[~(np.isnan(toxo_error))]  # Remove NaN values
    print(f"Global PyToxo error mean: {np.mean(pytoxo_errors)}")
    print(f"Global Toxo error mean: {np.mean(toxo_errors)}")
    print(f"Global PyToxo error sum: {np.sum(pytoxo_errors)}")
    print(f"Global Toxo error sum: {np.sum(toxo_errors)}")
    print(
        f"Global error mean comparison (`Toxo/PyToxo`): {np.mean(toxo_errors) / np.mean(pytoxo_errors)}"
    )
    print(
        f"Global error sums comparison (`Toxo/PyToxo`): {np.sum(toxo_errors) / np.sum(pytoxo_errors)}"
    )

    print()  # Newline

    # Calculate global accumulated times
    pytoxo_times = merged["Time (s) avg. 3_x"]
    toxo_times = merged["Time (s) avg. 3_y"]
    print(f"Global PyToxo time mean (s): {np.mean(pytoxo_times)}")
    print(f"Global Toxo time mean (s): {np.mean(toxo_times)}")
    print(f"Global PyToxo time sum (s): {np.sum(pytoxo_times)}")
    print(f"Global Toxo time sum (s): {np.sum(toxo_times)}")
    print(
        f"Global time mean comparison (`Toxo/PyToxo`): {np.mean(toxo_times) / np.mean(pytoxo_times)}"
    )
    print(
        f"Global time sums comparison (`Toxo/PyToxo`): {np.sum(toxo_times) / np.sum(pytoxo_times)}"
    )

    # Plot a graph since original data frames. First add new label to distinguish both
    pytoxo_report["Tool"] = ["PyToxo"] * pytoxo_report.shape[0]
    toxo_report["Tool"] = ["Toxo"] * toxo_report.shape[0]
    concat = pd.concat([pytoxo_report, toxo_report])
    graph = sns.relplot(
        data=concat,
        y="Error",
        x="Order",
        col="Model",
        hue="Tool",
        # size="Error",
        kind="scatter",
        alpha=0.5,
    )
    plt.show()

    # Save the generated graphic
    if SAVE:
        output_graph_path = os.path.join(
            "test",
            "accuracy",
            "reports",
            f"{script_name}_{now}_graph.png",
        )
        fig = graph.figure
        fig.savefig(output_graph_path)
