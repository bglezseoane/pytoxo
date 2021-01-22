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

import numpy as np
from sympy import *

from pytoxo.genotype_probabilities import genotype_probabilities
from pytoxo.model import Model


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
        self._order = model.order
        self._variables = {
            str(model.variables[0]): values[0],
            str(model.variables[1]): values[1],
        }
        self._penetrance_values = substitution(
            model.penetrances, model.variables, values
        )

    def compute_prevalence(self, mafs: [float], gp: [float] = None) -> float:
        """Compute the prevalence of the penetrance table.

        Parameters
        ----------
        mafs : [float]
            Minor allele frequencies array.
        gp : [float], optional
            Genotype probabilities array.

        Returns
        -------
        float
            Prevalence of the penetrance table.
        """
        if not gp:
            gp = genotype_probabilities(mafs)

        return float(np.sum(np.multiply(self._penetrance_values, gp)))

    def compute_heritability(self, mafs: [float]) -> float:
        """Compute the heritability of the penetrance table.

        Parameters
        ----------
        mafs : [float]
            Minor allele frequencies array.

        Returns
        -------
        float
            Heritability of the penetrance table.
        """
        gp = genotype_probabilities(mafs)
        p = self.compute_prevalence(mafs, gp)
        return float(
            np.sum(
                np.multiply(np.power(np.subtract(self._penetrance_values, p), 2)),
                gp,
            )
            / (p * (1 - p))
        )
