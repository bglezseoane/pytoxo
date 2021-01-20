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

from pytoxo.model import Model
from sympy import *


class PTable:
    """Representation of a penetrance table."""

    _order = None  # Number of locus defined by the penetrance table
    _variables = []  # Values for the variables present in the original model
    _penetrance_values = []  # Array of symbolic penetrances values

    def __init__(self, model: Model, values: [int]):
        """Creates a penetrance table from a given PyToxo model and its
        variable values.

        Parameters
        ----------
        model : Model
            Model from which to create the penetrance table.
        values : [int]
             Value for each of the variables represented in model.
        """
        self._order = model._order
        self._variables = {
            str(model._variables[0]): values[0],
            str(model._variables[1]): values[1],
        }
        self._penetrance_values = substitution(
            model._penetraces, model._variables, values
        )
