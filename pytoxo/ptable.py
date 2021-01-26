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
import os

import sympy


class PTable:
    """Representation of a penetrance table."""

    _order = None  # Number of locus defined by the penetrance table
    _variables = []  # Values for the variables present in the original model
    _penetrance_values = []  # Array of symbolic penetrances values

    def __init__(
        self,
        model_order: int,
        model_variables: list[sympy.Symbol],
        model_penetrances: list[float],
        values: [int],
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
        model_penetrances : list[float]
            Penetrances from the PyToxo model from which to create the
            penetrance table.
        values : [int]
             Value for each of the variables represented in model.
        """
        self._order = model_order
        self._variables = {
            str(model_variables[0]): values[0],
            str(model_variables[1]): values[1],
        }
        self._penetrance_values = sympy.substitution(
            model_penetrances, model_variables, values
        )

    def write_to_file(
        self, filename: str, overwrite: bool = False, format: str = "csv"
    ) -> None:
        """Write the penetrance table into a file.

        Currently only CSV format is supported.

        Parameters
        ----------
        filename : str
            The full file name where writes the table, without extension.
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
        # Calculate extension
        if format == "csv":
            extension = ".csv"
        else:
            ValueError(f"Unsupported '{format}' format")
        # Calculate final filename
        filename = os.path.normpath(filename + extension)
        # Check final file name
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(filename)
        if os.path.exists(filename) and not os.path.isfile(filename):
            raise IsADirectoryError(filename)

        with open(filename) as f:
            raise NotImplementedError  # TODO: Implement
