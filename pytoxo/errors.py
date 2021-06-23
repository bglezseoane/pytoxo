# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python tool to calculate penetrance tables for 
# high-order epistasis models
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

"""PyToxo error definitions."""

from typing import Dict, Tuple, Union


class BadFormedModelError(Exception):
    """Representation of a parse error with an epistatic model."""

    def __init__(
        self,
        error_object: Union[str, Dict[str, str], Tuple[list, list]],
        cause: str = "Bad formed model",
    ):
        """Creates an exception of a parse error with an epistatic model.

        Parameters
        ----------
        error_object : Union[str, dict[str, str], (list, list)]
            The input model which causes the error in one of its possible forms.
        cause : str, optional
            The error cause during parsing process (default "Bad formed file").
        """
        self.error_object = error_object
        self.cause = cause
        self.message = f"Bad formed model with '{error_object}' due to: '{cause}'."


class ResolutionError(Exception):
    """Representation of a resolution error of PyToxo."""

    def __init__(
        self,
        cause: str = "PyToxo can not solve this model",
        model_name: str = None,
        equation: str = None,
        message: str = None,
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
        message : str, optional
            The message to print, directly passed here. Other arguments would be
            ignored if this one is used.
        """
        msg_header = "Resolution error: "
        self.model_name = model_name
        self.cause = cause
        self.equation = equation
        self.message = message

        if not self.message:
            if not self.model_name:
                model_msg_text = "Model"
            else:
                model_msg_text = f"Model '{self.model_name}'"
            if self.equation:
                self.message = (
                    f"{msg_header} {model_msg_text} due to '{cause}'. "
                    f"Problematic equation: '{self.equation}'."
                )
            else:
                self.message = f"{msg_header} '{cause}'."
        else:
            self.message = f"{msg_header} {self.message}."


class UnsolvableModelError(Exception):
    """Representation of an unsolvable model error of PyToxo."""

    def __init__(
        self, model_name: str = None, equation: str = None, message: str = None
    ):
        """Creates an exception of an unsolvable model.

        Parameters
        ----------
        model_name : str, optional
            The model given name.
        equation : str, optional
            The equation that cause the resolution error, as string. Default is
            none.
        message : str, optional
            The message to print, directly passed here. Other arguments would be
            ignored if this one is used.
        """
        msg_header = "Unsolvable model:"
        self.model_name = model_name
        self.equation = equation
        self.message = message

        if not self.message:
            if self.model_name:
                model_msg_text = "model"
            else:
                model_msg_text = f"model ('{self.model_name}')"
            if self.equation:
                self.message = (
                    f"{msg_header} This {model_msg_text} has not solution. "
                    f"Equation: '{self.equation}'."
                )
            else:
                self.message = f"{msg_header} This {model_msg_text} has not solution."
        else:
            self.message = f"{msg_header} {self.message}."


class GenericCalculationError(Exception):
    """Representation of a generic calculation error of PyToxo."""

    def __init__(self, function: str = None):
        """Creates a generic calculation error of PyToxo.

        Parameters
        ----------
        function : str, optional
            The name of the function when the error occurs.
        """
        msg_skeleton = "Calculation error produced in: '{}'"
        self.function = function

        if not self.function:
            self.message = "Calculation error produced"
        else:
            self.message = msg_skeleton.format(self.function)


class GUIUnsupportedPlatformError(Exception):
    """Representation of a GUI platform support error of PyToxo."""

    def __init__(self, platform: str):
        """Creates a GUI platform support error of PyToxo.

        Parameters
        ----------
        platform : str
            The name of the platform where the GUI is used: Darwin, Linux,
            Windows or other unsupported one.
        """
        if platform == "Linux":
            self.message = (
                f"Your plaform (detected {platform}) supports PyToxo's GUI, "
                f"but something is wrong... Please, check you have Tk and "
                f"its Python binders installed and correctly reachable. "
                f"You can usually install them with the command:"
                f"\n\n\tsudo apt install python3-tk"
            )
        elif platform == "Darwin" or platform == "Windows":
            self.message = (
                f"Your plaform (detected {platform}) supports PyToxo's GUI, "
                f"but something is wrong... Please, check you have Tk and "
                f"its Python binders installed and correctly reachable."
            )
        else:
            self.message = (
                f"Your platform (detected {platform}) does not support the "
                f"PyToxo's GUI. Please, try the CLI instead."
            )
