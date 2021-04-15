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

"""Command line interface entry point."""

import argparse
import sys

import termcolor

import pytoxo
import pytoxo.errors

# Some predefined printing chunks
ok_hd = f"{termcolor.colored('[OK]', 'green', attrs=['bold'])}:"
warning_hd = f"{termcolor.colored('[WARNING]', 'yellow', attrs=['bold'])}:"
error_hd = f"{termcolor.colored('[ERROR]', 'red', attrs=['bold'])}:"


def main():
    # Define the parser
    parser = argparse.ArgumentParser(
        prog="pytoxo",
        description="PyToxo: A Python library for calculating penetrance "
        "tables of any bivariate epistasis model.",
        epilog="Copyright 2021 Borja González Seoane.",
        exit_on_error=False,
    )
    parser.add_argument(
        "model",
        metavar="MODEL",
        type=str,
        nargs=1,
        help="The path to the epistatic model CSV file.",
    )
    parser.add_argument(
        "maf",
        metavar="MAF",
        type=float,
        nargs=1,
        help="The MAF to use.",
    )
    prev_or_her_flag = parser.add_mutually_exclusive_group(required=True)
    prev_or_her_flag.add_argument(
        "--max_prev",
        action="store_true",
        help="This flag set PyToxo to maximize the prevalence, so the "
        "'PREV_OR_HER' argument will be considered as heritability.",
    )
    prev_or_her_flag.add_argument(
        "--max_her",
        action="store_true",
        help="This flag set PyToxo to maximize the prevalence, so the "
        "'PREV_OR_HER' argument will be considered as prevalence.",
    )
    parser.add_argument(
        "prev_or_her_arg",
        metavar="PREV_OR_HER",
        type=float,
        nargs=1,
        help="The heritability or prevalence to use, depending of the used flag.",
    )

    # Try to parse
    try:
        args = parser.parse_args()
    except argparse.ArgumentError:
        print(f"{error_hd} Argument parsing error.")
        sys.exit(1)

    # Get arguments
    model = args.model[0]
    maf = args.maf[0]
    if args.max_prev:
        calc_target = pytoxo.Model.find_max_prevalence_table
    else:
        calc_target = pytoxo.Model.find_max_heritability_table
    prev_or_her = args.prev_or_her_arg[0]

    # Try to build model
    try:
        model = pytoxo.Model(filename=model)
    except IOError as e:
        print(f"{error_hd} {e.message}")
        sys.exit(1)
    except pytoxo.errors.BadFormedModelError as e:
        print(f"{error_hd} {e.message}")
        sys.exit(1)

    # Try to solve the model
    try:
        ptable = calc_target(model, [maf] * model.order, prev_or_her)
    except pytoxo.errors.ResolutionError as e:
        print(f"{error_hd} {e.message}")
        sys.exit(1)
    except pytoxo.errors.UnsolvableModelError as e:
        print(f"{error_hd} {e.message}")
        sys.exit(1)

    # Print the penetrance table to can integrate it in a typical shell routine
    ptable.print_table()


if __name__ == "__main__":
    main()
