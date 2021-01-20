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
from sympy import Symbol, sympify, SympifyError

from pytoxo.errors import ModelCSVParsingError


class Model:
    """Representation of an epistasis model."""

    _name = None  # Name of the model
    _order = None  # Number of locus defined in the model
    _penetrances = []  # Array of symbolic expressions representing the epistatic model
    _variables = []  # List of symbolic variables used throughout the model

    def __init__(self, filename: str) -> None:
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
        ModelCSVParsingError
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
                raise ModelCSVParsingError(filename, "File without content.")

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
                    raise ModelCSVParsingError(
                        filename, "All the genotypes should have the same order."
                    )
            # Assert first members length is odd
            if len_fst_member % 2 != 0:
                raise ModelCSVParsingError(filename, "Bad genotype specification.")
            # Save the order
            self._order = len_fst_member // 2

            # Save the penetrances as symbolic expressions
            try:
                self._penetrances = [sympify(snd_member) for snd_member in snd_members]
            except SympifyError:
                raise ModelCSVParsingError(filename, "Bad penetrance specification.")

            # Save the variables of the used expressions
            for snd_member in snd_members:
                self._variables.append([Symbol(i) for i in snd_member if i.isalpha()])
        except IOError as e:
            raise e
        except ModelCSVParsingError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise ModelCSVParsingError(filename)
