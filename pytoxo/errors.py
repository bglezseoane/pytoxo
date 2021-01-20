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

"""PyToxo error definitions."""


class ModelCSVParsingError(Exception):
    """Representation of a parsing error with a model CSV file."""

    def __init__(self, filename: str, cause: str = "Bad formed file"):
        """Creates an exception of a parsing error with a model CSV file.

        Parameters
        ----------
        filename : str
            The path of the text file with the model CSV specification.
        cause : str
            The error cause during parsing process (default "Bad formed file").
        """
        self.filename = filename
        self.cause = cause
        self.message = f"ModelCSVParsingError with file '{filename}' due to: '{cause}'"
        super().__init__(self.message)
