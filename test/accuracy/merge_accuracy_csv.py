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

import pandas as pd  # TODO: Add to requirements

# ####################### EDIT HERE #######################
PYTOXO_REPORT_FILENAME = "test_variables_accuracy_2022-02-08_03-31-41.csv"
TOXO_REPORT_FILENAME = "build_toxo_accuracy_contrast_2022-02-08_17-51-47.csv"
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
