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

import numpy
import sympy

import pytoxo.calculations
import pytoxo.errors


class PTable:
    """Representation of a penetrance table."""

    def __init__(
        self,
        model_order: int,
        model_genotypes: list[str],
        model_penetrances: list[sympy.Expr],
        values: dict[sympy.Symbol, float],
        mafs: list[float] = None,
        model_name: str = None,
    ):
        """Creates a penetrance table from a given PyToxo model defined by
        its variables and penetrances, and its variable values.

        Parameters
        ----------
        model_order : int
            Order from the PyToxo model from which to create the penetrance
            table.
        model_genotypes : list[str]
            List of the genotype definitions like generated by
            `_calculate_genotypes` PyToxo model method.
        model_penetrances : list[sympy.Expr]
            Penetrances from the PyToxo model from which to create the
            penetrance table.
        values : dict[sympy.Symbol, float]
             Value for each of the variables represented in model, typically
             `x` and `y`. Keys are the own variables as `sympy.Symbol`.
        mafs : list[float], optional
             The MAFs used in the generation of the table. Optional
             parameter, only necessary to save then the table using GAMETES
             format.
        model_name: str, optional
            Optional value of the model which compound this table,
            to identify it easily.
        """
        self._order = model_order
        self._genotypes = model_genotypes
        self._penetrance_values = [
            p.subs(values) for p in model_penetrances
        ]  # Try to substitute `y` in expression `x` does not cause errors, simply are ignored
        self._values = list(values.values())  # Only used to save to GAMETES
        self._model_name = model_name
        self._mafs = mafs  # Only used to save to GAMETES
        self._prevalence = None  # Only used to save to GAMETES
        self._heritability = None  # Only used to save to GAMETES

    ########################################
    # Getters and setters for properties

    def _get_model_name(self) -> str:
        return self._model_name

    def _set_model_name(self, model_name: str) -> None:
        self._model_name = model_name

    model_name = property(_get_model_name, _set_model_name)

    def _get_mafs(self) -> list[float]:
        return self._mafs

    def _set_mafs(self, mafs: list[float]) -> None:
        self._mafs = mafs

    mafs = property(_get_mafs, _set_mafs)

    def _get_prevalence(self) -> float:
        return self._prevalence

    def _set_prevalence(self, prevalence: float) -> None:
        self._prevalence = prevalence

    prevalence = property(_get_prevalence, _set_prevalence)

    def _get_heritability(self) -> float:
        return self._heritability

    def _set_heritability(self, heritability: float) -> None:
        self._heritability = heritability

    heritability = property(_get_heritability, _set_heritability)

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
            hash(self._order)
            + hash(str(self._genotypes))
            + hash(str(self._penetrance_values))
            + hash(self._values)
            + hash(self._model_name)
            + hash(self._mafs)
            + self._prevalence
            + self._heritability
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

    def _compound_table_as_gametes(self) -> str:
        """Compound the penetrance table as GAMETES format, to save it into a
        file.

        Attributes `mafs`, `prevalence` and `heritability` must be filled.

        Returns
        -------
        str
            The penetrance table formatted as a text string.

        Raises
        ------
        GenericCalculationError
            Failing to calculate prevalence or heritability to save as GAMETES.
        """
        """Calculate prevalence and heritability. Do here to optimize,
            because they are only needed to save as GAMETES format"""
        try:
            self._prevalence = pytoxo.calculations.compute_prevalence(
                penetrances=self._penetrance_values,
                mafs=self._mafs,
                model_order=self._order,
            )
            self._heritability = pytoxo.calculations.compute_heritability(
                penetrances=self._penetrance_values,
                mafs=self._mafs,
                model_order=self._order,
            )
        except pytoxo.errors.GenericCalculationError as e:
            raise pytoxo.errors.GenericCalculationError(
                "It is no possible to calculate the "
                f"prevalence or the heritability using {e.function}."
                "Are you using the original used MAFs?"
            )

        gametes_skeleton = (
            "Attribute names:\t{}\nMinor allele "
            "frequencies:\t{}\nx: {}\ny: {}\nPrevalence: {"
            "}\nHeritability: {}\n\nTable:\n\n{} "
        )

        # Prepare fields to fill
        attribute_names = "\t".join([f"P{i}" for i in range(0, self.order)])
        mafs = "\t".join([f"{'{:.3f}'.format(maf)}" for maf in self._mafs])
        x = str(self._values[0])
        y = str(self._values[1])
        prev = str(self._prevalence.evalf())
        her = str(self._heritability.evalf())
        table = ", ".join([str(pen) for pen in self._penetrance_values])

        # Generate lines of the file with genotypes and its penetrances
        return gametes_skeleton.format(attribute_names, mafs, x, y, prev, her, table)

    def print_table(self) -> None:
        """Print the penetrance as raw text."""
        print(self._compound_table_as_text())

    def write_to_file(
        self, filename: str, overwrite: bool = False, format: str = "csv"
    ) -> None:
        """Write the penetrance table into a file.

        Currently formats CSV and GAMETES are available. To use GAMETES,
        this object has to be filled the attributes `mafs`, `prevalence`
        and `heritability`.

        Parameters
        ----------
        filename : str
            The full file name where writes the table.
        overwrite : bool (default False)
            A flag that should be passed as true to overwrite the final file
            if it already exists.
        format : str
            The format of the final file. Currently CSV format (`csv`
            flag) and GAMETES (`gametes` flag) are supported. Default is CSV.

        Raises
        ------
        ValueError
            If unsupported format tentative.
        FileExistsError
            If the `filename` file already exists and the `overwrite` has not be
            passed as true.
        IsADirectoryError
            If `filename` is a existent directory.
        GenericCalculationError
            Failing to calculate prevalence or heritability to save as GAMETES.
        """
        supported_formats = ["csv", "gametes"]

        # Input handling and checks
        format = format.lower()  # Defer
        if format not in supported_formats:
            ValueError(f"Unsupported '{format}' format")

        # Check if is possible to use GAMETES
        if format == "gametes" and (
            not self._mafs or not self._prevalence or not self._heritability
        ):
            raise ValueError(
                f"The '{format}' format requires to be filled the "
                f"attributes `mafs`, `prevalence`and `heritability`."
            )

        # Calculate final filename
        filename = os.path.normpath(filename)
        # Check final file name
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(filename)
        if os.path.exists(filename) and not os.path.isfile(filename):
            raise IsADirectoryError(filename)

        # Generate the table
        if format == "csv":
            table = self._compound_table_as_text()
        else:
            table = self._compound_table_as_gametes()

        # Write file
        if overwrite:
            os.remove(filename)
        with open(filename, "x") as f:
            f.write(table)
