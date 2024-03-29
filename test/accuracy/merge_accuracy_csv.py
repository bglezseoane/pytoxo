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

import numpy as np
import pandas as pd  # TODO: Add to requirements

# ####################### EDIT HERE #######################
PYTOXO_REPORT_FILENAME = (
    "test_variables_accuracy_max_p_2022-02-16_12-30-11.csv"  # Max. p.
)
# PYTOXO_REPORT_FILENAME = (
#     "test_variables_accuracy_max_h_2022-02-16_17-08-02.csv"  # Max. h.
# )
TOXO_REPORT_FILENAME = (
    "build_toxo_accuracy_contrast_max_p_2022-02-16_15-54-18.csv"  # Max. p.
)
# TOXO_REPORT_FILENAME = (
#     "build_toxo_accuracy_contrast_max_h_2022-02-16_18-11-42.csv"  # Max. h.
# )

SAVE = True
SOLUBILITY_CONTRIBUTION = True

# Comment or uncomment firsts or seconds of each pair
PREV_OR_HER_STR = "Prevalence"
# PREV_OR_HER_STR = "Heritability"
PREV_OR_HER_LETTER = "p"
# PREV_OR_HER_LETTER = "h"

MODELS_SCOPE = [
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
    "multiplicative_5",
    "threshold_2",
    "threshold_3",
    "threshold_4",
    "threshold_5",
    "threshold_6",
    "threshold_7",
    "threshold_8",
]  # Other models are ignored
# #########################################################

# Calculate some necessary opposites
if PREV_OR_HER_STR == "Prevalence":
    prev_or_her_str_op = "Heritability"
else:
    prev_or_her_str_op = "Prevalence"
if PREV_OR_HER_LETTER == "p":
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
    f"{script_name}_max_{PREV_OR_HER_LETTER}_{now}.csv",
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
    merged_left = pytoxo_report.merge(
        toxo_report,
        on=["Model", "Order", "MAF", f"{prev_or_her_str_op}"],
        how="left",
        indicator=True,
    )
    merged_right = pytoxo_report.merge(
        toxo_report,
        on=["Model", "Order", "MAF", f"{prev_or_her_str_op}"],
        how="right",
        indicator=True,
    )

    # Clean models out of scope
    models_toxo_df = merged_right.iloc[:, :2].drop_duplicates()
    # Different model-order combs. in Toxo data
    models_toxo = []
    for _, row in models_toxo_df.iterrows():
        models_toxo.append(f"{row['Model'].lower()}_{row['Order']}")
    models_out_of_scope = [
        [m.split("_")[0].capitalize(), m.split("_")[1]]
        for m in models_toxo
        if m not in MODELS_SCOPE
    ]
    # Remove from original data frame
    for model, order in models_out_of_scope:
        merged_right.drop(
            merged_right[
                (merged_right.Model == model) & (merged_right.Order == int(order))
            ].index,
            inplace=True,
        )

    # Achieve masks: only-Pytoxo, only_Toxo, both
    only_pytoxo_mask = merged_left["_merge"] == "left_only"
    both_mask = merged_left["_merge"] == "both"
    only_toxo_mask = merged_right["_merge"] == "right_only"

    # Prepare solubility contribution aux. data frames
    only_pytoxo = merged_left[only_pytoxo_mask].drop(
        "_merge", 1
    )  # Drop useless columns
    both = merged_left[both_mask].drop("_merge", 1)  # Drop useless columns
    only_toxo = merged_right[only_toxo_mask].drop("_merge", 1)  # Drop useless columns

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
                f"{script_name}_max_{PREV_OR_HER_LETTER}_only_pytoxo_{now}.csv",
            ),
            os.path.join(
                contribution_report_base_path,
                f"{script_name}_max_{PREV_OR_HER_LETTER}_both_{now}.csv",
            ),
            os.path.join(
                contribution_report_base_path,
                f"{script_name}_max_{PREV_OR_HER_LETTER}_only_toxo_{now}.csv",
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
pytoxo_times = merged["Time (s) avg. 3"]
toxo_times = merged["Time (s) avg."]
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
