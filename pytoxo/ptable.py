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

import numpy
import sympy


class PTable:
    """Representation of a penetrance table."""

    def __init__(
        self,
        model_order: int,
        model_penetrances: list[sympy.Expr],
        values: dict[sympy.Symbol, float],
        model_name: str = None,
    ):
        """Creates a penetrance table from a given PyToxo model defined by
        its variables and penetrances, and its variable values.

        Parameters
        ----------
        model_order : int
            Order from the PyToxo model from which to create the penetrance
            table.
        model_penetrances : list[sympy.Expr]
            Penetrances from the PyToxo model from which to create the
            penetrance table.
        values : dict[sympy.Symbol: float]
             Value for each of the variables represented in model, typically
             `x` and `y`. Keys are the own variables as `sympy.Symbol`.
        model_name: str, optional
            Optional value of the model which compound this table,
            to identify it easily.
        """
        self._order = model_order
        self._calculate_genotypes()
        self._penetrance_values = [
            p.subs(values) for p in model_penetrances
        ]  # Try to substitute `y` in expression `x` does not cause errors, simply are ignored
        self._model_name = model_name

    def _calculate_genotypes(self) -> None:
        """Calculate the list of genotypes to the given model attending to
        the model order.

        Uses default alphabetical sort like in the `Model` constructors.
        Capital letters first.
        """
        # Deduce alleles attending to the table order
        # TODO: Is `len(ascii_lowercase) < self._order` possible?
        letters = list(string.ascii_lowercase[: self._order])
        alleles = []
        for letter in letters:
            lower = letter
            upper = letter.upper()
            alleles.append([f"{upper}{upper}", f"{upper}{lower}", f"{lower}{lower}"])
        # Generate genotypes tracing the alleles cartesian product
        self._genotypes = [str("".join(p)) for p in list(itertools.product(*alleles))]

    ########################################
    # Getters and setters for properties

    def _get_model_name(self) -> str:
        return self._model_name

    def _set_model_name(self, model_name: str) -> None:
        self._model_name = model_name

    model_name = property(_get_model_name, _set_model_name)

    @property
    def order(self) -> int:
        return self._order

    @property
    def genotypes(self) -> list[str]:
        return self._genotypes

    @property
    def genotypes_as_numpy(self) -> numpy.array:
        return numpy.array(self._genotypes)

    @property
    def penetrance_values(self) -> list[float]:
        return self._penetrance_values

    @property
    def penetrance_values_as_numpy(self) -> numpy.array:
        return numpy.array(self._penetrance_values)

    ########################################

    def __hash__(self):
        return hash(
            hash(self._model_name)
            + hash(self._order)
            + hash(str(self._genotypes))
            + hash(str(self._penetrance_values))
        )

    def __eq__(self, other):
        return hash(self) == hash(other)

    def _compound_table_as_text(self) -> str:
        """Compound the penetrance table as text, to print it ot save into a
        file it.

        Returns
        -------
        str
            The penetrance table formatted as a text string.
        """
        # Generate lines of the file with genotypes and its penetrances
        lines = ""
        for genotype, penetrance in zip(self._genotypes, self._penetrance_values):
            lines += f"{''.join(genotype)},{penetrance}\n"

        return lines

    def print_table(self) -> None:
        """Print the penetrance as raw text."""
        print(self._compound_table_as_text())

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
        if not format == "csv":
            ValueError(f"Unsupported '{format}' format")
        # Calculate final filename
        filename = os.path.normpath(filename)
        # Check final file name
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(filename)
        if os.path.exists(filename) and not os.path.isfile(filename):
            raise IsADirectoryError(filename)

        # Write file
        with open(filename, "x") as f:
            f.write(self._compound_table_as_text())
