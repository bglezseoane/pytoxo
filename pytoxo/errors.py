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


class ModelCSVParseError(Exception):
    """Representation of a parse error with a model CSV file."""

    def __init__(self, filename: str, cause: str = "Bad formed file"):
        """Creates an exception of a parse error with a model CSV file.

        Parameters
        ----------
        filename : str
            The path of the text file with the model CSV specification.
        cause : str, optional
            The error cause during parsing process (default "Bad formed file").
        """
        self.filename = filename
        self.cause = cause
        self.message = f"ModelCSVParseError with file '{filename}' due to: '{cause}'"
        super().__init__(self.message)


class ResolutionError(Exception):
    """Representation of a resolution error of PyToxo."""

    def __init__(
        self,
        cause: str = "PyToxo can not solve this model",
        model_name: str = None,
        equation: str = None,
    ):
        """Creates an exception of a solving error with a PyToxo model.

        Parameters
        ----------
        cause : str, optional
            The cause during solve process (default "PyToxo can not solve
            this model").
        model_name : str, optional
            The model given name.
        equation : str, optional
            The equation that cause the resolution error, as string. Default is
            none.
        """
        self.model_name = model_name
        self.cause = cause
        self.equation = equation
        if not self.model_name:
            model_msg_text = "model"
        else:
            model_msg_text = f"model ('{self.model_name}')"
        if self.equation:
            self.message = (
                f"ResolutionError with {model_msg_text}: '{cause}'. "
                f"Problematic equation: '{self.equation}'."
            )
        else:
            self.message = f"ResolutionError: '{cause}'."
        super().__init__(self.message)


class UnsolvableModelError(Exception):
    """Representation of an unsolvable model error of PyToxo."""

    def __init__(self, model_name: str = None, equation: str = None):
        """Creates an exception of an unsolvable model.

        Parameters
        ----------
        model_name : str, optional
            The model given name.
        equation : str, optional
            The equation that cause the resolution error, as string. Default is
            none.
        """
        self.model_name = model_name
        self.equation = equation
        if self.model_name:
            model_msg_text = "model"
        else:
            model_msg_text = f"model ('{self.model_name}')"
        if self.equation:
            self.message = (
                f"UnsolvableModelError: This {model_msg_text} has not solution. "
                f"Equation: '{self.equation}'."
            )
        else:
            self.message = (
                f"UnsolvableModelError: This {model_msg_text} has not solution."
            )
        super().__init__(self.message)
