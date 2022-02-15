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

Enabling the flag `SOLUBILITY_CONTRIBUTION` this test also prepares a
contribution report for the solubility collection, taking advantage of this
data. This is only interesting handling reports (`PYTOXO_REPORT_FILENAME` and
`TOXO_REPORT_FILENAME`) with a big amount of cases within.
"""

import datetime
import os
import sys

import pandas as pd  # TODO: Add to requirements

# ####################### EDIT HERE #######################
PYTOXO_REPORT_FILENAME = "test_variables_accuracy_2022-02-08_03-31-41.csv"
TOXO_REPORT_FILENAME = "build_toxo_accuracy_contrast_2022-02-08_17-51-47.csv"
SAVE = True
SOLUBILITY_CONTRIBUTION = True

# Comment or uncomment firsts or seconds of each pair
prev_or_her_str = "Prevalence"
# prev_or_her_str = "Heritability"
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
merged = pytoxo_report.merge(
    toxo_report, on=["Model", "Order", "MAF", f"{prev_or_her_str_op}"]
)

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
    f"{script_name}_max_{prev_or_her_letter}_{now}.csv",
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
                    f"{prev_or_her_str_op}",
                    "Error PyToxo",
                    "Time (s) avg. 3 PyToxo",
                    "Error Toxo",
                    "Time (s) avg. Toxo\n",
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


if SOLUBILITY_CONTRIBUTION:
    # Merge reports again, but preserving cases only solvable by Pytoxo or by
    # Toxo
    merged = pytoxo_report.merge(
        toxo_report,
        on=["Model", "Order", "MAF", f"{prev_or_her_str_op}"],
        how="left",
        indicator=True,
    )

    # Achieve masks: only-Pytoxo, only_Toxo, both
    only_pytoxo_mask = merged["_merge"] == "left_only"
    both_mask = merged["_merge"] == "both"
    only_toxo_mask = merged["_merge"] == "right_only"

    # Prepare solubility contribution aux. data frames
    only_pytoxo = merged[only_pytoxo_mask].drop("_merge", 1)  # Drop useless columns
    both = merged[both_mask].drop("_merge", 1)  # Drop useless columns
    only_toxo = merged[only_toxo_mask].drop("_merge", 1)  # Drop useless columns

    if SAVE:
        # Prepare output paths
        contribution_report_base_path = os.path.join(
            "test",
            "solubility",
            "reports",
        )
        contribution_report_paths = [
            os.path.join(
                contribution_report_base_path,
                f"{script_name}_max_{prev_or_her_letter}_only_pytoxo_{now}.csv",
            ),
            os.path.join(
                contribution_report_base_path,
                f"{script_name}_max_{prev_or_her_letter}_both_{now}.csv",
            ),
            os.path.join(
                contribution_report_base_path,
                f"{script_name}_max_{prev_or_her_letter}_only_toxo_{now}.csv",
            ),
        ]

        # Save contribution data frames as CSV
        for filename, dataframe in zip(
            contribution_report_paths, [only_pytoxo, both, only_toxo]
        ):
            with open(filename, "w") as f:  # Patch to set custom headers
                f.write(
                    ";".join(
                        [
                            "Model",
                            "Order",
                            "MAF",
                            f"{prev_or_her_str_op}",
                            "Error PyToxo",
                            "Time (s) avg. 3 PyToxo",
                            "Error Toxo",
                            "Time (s) avg. Toxo\n",
                        ]
                    )
                )
            dataframe.to_csv(
                filename,
                header=False,
                index=False,
                sep=";",
                mode="a",
            )
