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
import statistics
import typing

import sympy
import timeout_decorator

import pytoxo.errors
import pytoxo.ptable
import pytoxo.calculations


class Model:
    """Representation of an epistasis model."""

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
        self._name = None  # Name of the model
        self._order = None  # Number of locus defined in the model
        self._penetrances = (
            []
        )  # Array of symbolic expressions representing the epistatic model
        self._variables = []  # List of symbolic variables used throughout the model

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
            fst_members = [line.split(",")[0].strip() for line in lines]
            snd_members = [line.split(",")[1].strip() for line in lines]

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
            all_variables = []
            for snd_member in snd_members:
                all_variables.append(
                    sympy.symbols([i for i in snd_member if i.isalpha()])
                )

            # Reduce the variables to avoid replication
            all_variables = [
                item for subl in all_variables for item in subl
            ]  # Flat list
            different_variables = []
            for variable in all_variables:
                """Not use other more direct method like
                `list(set(all_variables))` because the variable apparition
                order must be preserved"""
                if variable not in different_variables:
                    different_variables.append(variable)
                    # TODO: Test with other variable names to assert...

            # Check support: only 2 different variables are supported by PyToxo
            if not 0 < len(different_variables) <= 2:
                raise pytoxo.errors.ModelCSVParsingError(
                    filename,
                    "Only models with 1 or 2 different variables are supported.",
                )
            else:
                self._variables = different_variables
        except IOError as e:
            raise e
        except pytoxo.errors.ModelCSVParsingError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise pytoxo.errors.ModelCSVParsingError(filename)

    ########################################
    # Getters and setters

    def get_name(self) -> str:
        return self._name

    def set_name(self, name: str) -> None:
        self._name = name

    name = property(get_name, set_name)

    @property
    def order(self) -> int:
        return self._order

    @property
    def penetrances(self) -> list[typing.Union[float, sympy.Expr]]:
        return self._penetrances

    @property
    def variables(self) -> list[sympy.Symbol]:
        return self._variables

    ########################################

    def _max_penetrance(self) -> sympy.Expr:
        """Returns the largest polynomial from all penetrance expressions, for
        any real and positive value of the two variables."""

        unique_penetrances = list(set(self._penetrances))

        # Return the maximum penetrance expression with the largest degree
        # TODO: Temporal patch. Revise math approach to achieve an stable
        #  solution to this step...
        degrees = [max(sympy.Poly(p).degree_list()) for p in unique_penetrances]
        return unique_penetrances[degrees.index(max(degrees))]

    def _solve(self, constraints: list[sympy.Eq]) -> pytoxo.ptable.PTable:
        """Solves the equation system formed by the provided equations.

        Parameters
        ----------
        constraints : list[sympy.Eq]
            Input constraints that define the equation.

        Returns
        -------
        pytoxo.ptable.PTable
            The equation solution.
        """
        # TODO: Consider add assumptions to vars real and greather than 0
        """Sometimes the solver does not find the correct solution, so the 
        next lines calculate an initial solution and then a confirmation one. 
        If the two are equivalent, the solution is considered as the correct. 
        If not, more calls to the solver are required to achieve a 75 % of 
        votation for the final one."""

        sol_minimum_acceptance_frequency = 0.75
        max_steps = 5
        solve_timeout = 10
        risky = False

        def try_to_solve() -> list[tuple[float]]:
            """Tries solve the given constraints with Sympy and returns
            solutions, filtering unreal and negative ones."""

            @timeout_decorator.timeout_decorator.timeout(
                solve_timeout, timeout_exception=StopIteration
            )
            def solver_call():
                """Help function which encapsulate the call to the solver
                with a timeout decorator to abort futile executions."""
                return sympy.solve(
                    constraints,
                    self._variables[0],
                    self._variables[1],
                    manual=True,
                    rational=False,
                )

            # Try to solve the system within the setting timeout
            try:
                sols = solver_call()
            except StopIteration:
                return []

            # Discard unreal solutions
            sols = [s for s in sols if s[0].is_real and s[1].is_real]
            # Discard negative solutions
            return [s for s in sols if s[0] > 0 and s[1] > 0]

        def check_solutions(sols) -> tuple[bool, typing.Union[tuple[float], None]]:
            """Checks if there are a valid solution in the input solutions
            list (a valid solution is one that has been reached at least with
            the `sol_minimum_acceptance_frequency` frequency). Returns true of
            false about the above check and also the mode solution. Returns
            false and none if a empty list is passed. Returns false and the
            unique solution is a list with only one solution is passed."""
            if not sols:
                return False, None
            elif len(sols) == 1:
                return False, sols[0]
            else:
                mode_sol = statistics.mode(sols)
                if sols.count(mode_sol) >= sol_minimum_acceptance_frequency * len(sols):
                    return True, mode_sol
                else:
                    return False, mode_sol

        sols = []
        # First tentative to achieve the solution
        for _ in range(max_steps):
            sols.extend(try_to_solve())
            if sols:
                sol = sols[0]  # Simply use the first one, if exists
                break

        # If not in risky mode, confirm the solution
        if not risky:
            ok = False
            while not ok:
                sols.extend(try_to_solve())
                ok, sol = check_solutions(sols)
                if len(sols) >= max_steps:  # Protection against infinite loop
                    break

        # TODO: `solve` return check and raising
        if not sol:
            raise ValueError("PyToxo can not solve this model.")

        # Return the final achieved solution
        return pytoxo.ptable.PTable(
            model_order=self._order,
            model_penetrances=self._penetrances,
            values={self._variables[0]: sol[0], self._variables[1]: sol[1]},
        )

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
        c1 = pytoxo.calculations.compute_heritability(self._penetrances, mafs) - h
        c2 = self._max_penetrance() - sympy.Integer(1)
        return self._solve(constraints=[c1, c2])

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
        c1 = pytoxo.calculations.compute_prevalence(self._penetrances, mafs) - p
        c2 = self._max_penetrance() - sympy.Integer(1)
        return self._solve(constraints=[c1, c2])
