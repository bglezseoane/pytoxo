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

"""Epistasis model definition."""

import os

import numpy as np
import sympy

import pytoxo.errors
import pytoxo.ptable
import pytoxo.calculations


class Model:
    """Representation of an epistasis model."""

    _name = None  # Name of the model
    _order = None  # Number of locus defined in the model
    _penetrances = []  # Array of symbolic expressions representing the epistatic model
    _variables = []  # List of symbolic variables used throughout the model

    def __init__(self, filename: str):
        """Reads the model from its text representation in `file` and inits
        an object with its data.

        The input model must be formatted as a plain CSV, with each line of the
        file corresponding to a row of the model. The rows are made of the
        genotype definition and the probability associated with the given
        genotype, separated by a comma. Probability is expressed as a function
        of two variables, represented as alphabetical characters. Empty lines,
        as well as lines starting with '#' will be ignored.

        Parameters
        ----------
        filename : str
            The path of the text file with the model.

        Raises
        ------
        pytoxo.errors.ModelCSVParsingError
            If the parsing tentative fails due to the file is not well formed.
        IOError
           If the input file is not found or the reading tentative fails due
           to other unexpected operative system level cause.
        """
        try:
            with open(filename, "r") as f:
                lines = f.readlines()

            # Discard comments and empty lines
            lines = [line for line in lines if not line.startswith("#")]
            lines = [line for line in lines if line.strip()]

            # Check not empty file
            if not lines:
                raise pytoxo.errors.ModelCSVParsingError(
                    filename, "File without content."
                )

            # Split lines around the comma (',')
            fst_members = [line.split(",")[0] for line in lines]
            snd_members = [line.split(",")[1] for line in lines]

            # Save the name of the model
            self._name = os.path.basename(filename).split(".")[0]

            # Catch the order of the model
            len_fst_member = len(fst_members[0])
            # Assert all the first members have the same length
            for fst_member in fst_members:
                if len(fst_member) != len_fst_member:
                    raise pytoxo.errors.ModelCSVParsingError(
                        filename, "All the genotypes should have the same order."
                    )
            # Assert first members length is odd
            if len_fst_member % 2 != 0:
                raise pytoxo.errors.ModelCSVParsingError(
                    filename, "Bad genotype specification."
                )
            # Save the order
            self._order = len_fst_member // 2

            # Save the penetrances as symbolic expressions
            try:
                self._penetrances = [
                    sympy.sympify(snd_member) for snd_member in snd_members
                ]
            except sympy.SympifyError:
                raise pytoxo.errors.ModelCSVParsingError(
                    filename, "Bad penetrance specification."
                )

            # Save the variables of the used expressions
            for snd_member in snd_members:
                self._variables.append(
                    [sympy.Symbol(i) for i in snd_member if i.isalpha()]
                )
            # Check support: only 2 variables are supported by PyToxo
            for vars in self._variables:
                if not len(vars) <= 2:
                    raise pytoxo.errors.ModelCSVParsingError(
                        filename, "Only models with 2 variables are supported."
                    )
        except IOError as e:
            raise e
        except pytoxo.errors.ModelCSVParsingError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise pytoxo.errors.ModelCSVParsingError(filename)

    ########################################
    # Getters and setters

    def get_name(self) -> None:
        return self._name

    def set_name(self, name: str) -> None:
        self._name = name

    name = property(get_name, set_name)

    @property
    def order(self) -> None:
        return self._order

    @property
    def penetrances(self) -> None:
        return self._penetrances

    @property
    def variables(self) -> None:
        return self._variables

    ########################################

    def _max_penetrance(self) -> None:
        """Returns the largest polynomial from all penetrance expressions, for
        any real and positive value of the two variables."""
        p = np.transpose(
            np.unique(
                np.array(filter(lambda x: sympy.degree(x) > 0, self._penetrances))
            )
        )

        """Convert model variables to constraints to extend them and build the
        final expression"""
        var_constraints = []
        for var in self._variables:
            var_constraints.append(sympy.sympify(sympy.Symbol(var) >= 0))

        tmp_max = p[1]  # Temporary maximum
        for i in p[2:]:
            constraints = var_constraints  # Load computed vars constraints
            # TODO: Check only real numbers are returned in `solve`
            s_x, _ = sympy.solve(
                constraints.extend(
                    [
                        sympy.sympify(sympy.Symbol(p) >= 0),
                        sympy.sympify(sympy.Symbol(p) <= 1),
                        sympy.sympify(sympy.Symbol(tmp_max) > sympy.Symbol(i)),
                    ]
                )
            )
            if not s_x:
                tmp_max = i
        return tmp_max

    def _solve(
        self, constraints: list[sympy.core.relational.Relational]
    ) -> tuple[float]:
        """Solves the equation system formed by the provided equations.

        Parameters
        ----------
        constraints : list[sympy.core.relational.Relational]
            Input constraints that define the equation.

        Returns
        -------
        tuple[float]
            The equation solution.
        """
        # TODO: Consider add assumptions to vars real and greather than 0
        [s_x, s_y] = sympy.solve(constraints, self._variables)
        # TODO: `solve` return check and raising
        return s_x, s_y

    def find_max_prevalence(self, mafs: list[float], h: float) -> pytoxo.ptable.PTable:
        """Computes the table whose prevalence is maximum for the given MAFs
        and heritability, and returns it within a `PTable` object.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        h : float
            Heritability of thr table.

        Returns
        -------
        pytoxo.ptable.PTable
            Penetrance table obtained within a `PTable` object.
        """
        c1 = sympy.sympify(
            pytoxo.calculations.compute_heritability(self._penetrances, mafs) == h
        )
        c2 = sympy.sympify(self._max_penetrance() == 1)
        [s_x, s_y] = self._solve(constraints=[c1, c2])
        return pytoxo.ptable.PTable(
            model_order=self._order,
            model_variables=self._variables,
            model_penetrances=self._penetrances,
            values=[s_x, s_y],
        )

    def find_max_heritability(
        self, mafs: list[float], p: float
    ) -> pytoxo.ptable.PTable:
        """Computes the table whose heritability is maximum for the given MAFs
        and prevalence, and returns it within a `PTable` object.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        p : float
            Prevalence of the table.

        Returns
        -------
        pytoxo.ptable.PTable
            Penetrance table obtained within a `PTable` object.
        """
        c1 = sympy.sympify(
            pytoxo.calculations.compute_prevalence(self._penetrances, mafs) == p
        )
        c2 = sympy.sympify(self._max_penetrance() == 1)
        [s_x, s_y] = self._solve(constraints=[c1, c2])
        return pytoxo.ptable.PTable(
            model_order=self._order,
            model_variables=self._variables,
            model_penetrances=self._penetrances,
            values=[s_x, s_y],
        )
