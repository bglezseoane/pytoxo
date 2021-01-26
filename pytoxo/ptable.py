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
