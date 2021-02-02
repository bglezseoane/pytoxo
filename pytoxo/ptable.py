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

"""Penetrance table definition."""

import itertools
import os
import string

import sympy


class PTable:
    """Representation of a penetrance table."""

    def __init__(
        self,
        model_order: int,
        model_variables: list[sympy.Symbol],
        model_penetrances: list[sympy.Expr],
        values: list[float],
    ):
        """Creates a penetrance table from a given PyToxo model defined by
        its variables and penetrances, and its variable values.

        Parameters
        ----------
        model_order : int
            Order from the PyToxo model from which to create the penetrance
            table.
        model_variables : list[sympy.Symbol]
            Variables from the PyToxo model from which to create the penetrance
            table.
        model_penetrances : list[sympy.Expr]
            Penetrances from the PyToxo model from which to create the
            penetrance table.
        values : list[float]
             Value for each of the variables represented in model.
        """
        self._order = model_order
        self._variables = model_variables
        substitution = {(self._variables[0]): values[0], self._variables[1]: values[1]}
        self._penetrance_values = [
            p.subs(substitution) for p in model_penetrances
        ]  # Try to substitute `y` in expression `x` does not cause errors, simply are ignored

    def write_to_file(
        self, filename: str, overwrite: bool = False, format: str = "csv"
    ) -> None:
        """Write the penetrance table into a file.

        Currently only CSV format is supported.

        Parameters
        ----------
        filename : str
            The full file name where writes the table.
        overwrite : bool (default False)
            A flag that should be passed as true to overwrite the final file
            if it already exists.
        format : str
            The format of the final file. Currently only CSV format (`csv`
            flag) is supported.

        Raises
        ------
        ValueError
            If unsupported format tentative.
        FileExistsError
            If the `filename` file already exists and the `overwrite` has not be
            passed as true.
        IsADirectoryError
            If `filename` is a existent directory.
        """
        # Input handling and checks
        if format == "CSV":  # Deference fix
            format = "csv"
        else:
            ValueError(f"Unsupported '{format}' format")
        # Calculate final filename
        filename = os.path.normpath(filename)
        # Check final file name
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(filename)
        if os.path.exists(filename) and not os.path.isfile(filename):
            raise IsADirectoryError(filename)

        # Generate genotypes column with genotype names
        # Deduce alleles attending to the table order
        # TODO: Is `len(ascii_lowercase) < self._order` possible?
        letters = list(string.ascii_lowercase[: self._order])
        alleles = []
        for letter in letters:
            lower = letter
            upper = letter.upper()
            alleles.append([f"{upper}{upper}", f"{upper}{lower}", f"{lower}{lower}"])
        # Generate genotypes tracing the alleles cartesian product
        genotypes = list(itertools.product(*alleles))

        # Generate lines of the file with genotypes and its penetrances
        lines = []
        for genotype, penetrance in zip(genotypes, self._penetrance_values):
            lines.append(f"{''.join(genotype)},{penetrance}\n")

        # Write file
        with open(filename, "x") as f:
            f.writelines(lines)
