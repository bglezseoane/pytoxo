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
import typing

import mpmath
import sympy
import timeout_decorator

import pytoxo.calculations
import pytoxo.errors
import pytoxo.ptable

_MPMATH_DEFAULT_DPS = 15
_TOLERABLE_SOLUTION_ERROR_BASE_DELTA = 1e-16  # It is fitted then to model's order


class Model:
    """Representation of an epistasis model."""

    def __init__(self, filename: str):
        """Init the object.

        Uses the model from its text representation in `filename` file and
        inits an object with its data.

        The input model must be formatted as a plain CSV, with each line of the
        file corresponding to a row of the model. The rows are made of the
        genotype definition and the probability associated with the given
        genotype, separated by a comma.

        Probability is expressed as a function of two variables, represented
        as alphabetical characters. Only two variables are supported at most.
        Any pair of alphabetic characters can be used to represent these
        variables, but they must be the same two names for the entire model
        file.

        Empty lines, as well as lines starting with '#' will be ignored.

        Parameters
        ----------
        filename : str
            The path of the text file with the model.

        Raises
        ------
        pytoxo.errors.ModelCSVParseError
            If the parse tentative fails due to the file is not well formed.
        IOError
           If the input file is not found or the reading tentative fails due
           to other unexpected operative system level cause.
        """
        self._name = None  # Name of the model
        self._order = None  # Number of locus defined in the model
        self._penetrances = (
            []
        )  # List of symbolic expressions representing the epistatic model
        self._variables = []  # List of symbolic variables used throughout the model
        self._tolerable_solution_error_delta = None

        # Delegates parse to the `parse` method
        self._parse(filename)

    def _parse(self, filename: str) -> None:
        """Takes the responsibility of the initializer to parse the model file.

        Reads the model from its text representation in `filename` file and
        inits the object with its data.

        The input model must be formatted as a plain CSV, with each line of the
        file corresponding to a row of the model. The rows are made of the
        genotype definition and the probability associated with the given
        genotype, separated by a comma.

        Probability is expressed as a function of two variables, represented
        as alphabetical characters. Only two variables are supported at most.
        Any pair of alphabetic characters can be used to represent these
        variables, but they must be the same two names for the entire model
        file.

        Empty lines, as well as lines starting with '#' will be ignored.

        Parameters
        ----------
        filename : str
            The path of the text file with the model.

        Raises
        ------
        pytoxo.errors.ModelCSVParseError
            If the parse tentative fails due to the file is not well formed.
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
                raise pytoxo.errors.ModelCSVParseError(
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
                    raise pytoxo.errors.ModelCSVParseError(
                        filename, "All the genotypes should have the same order."
                    )
            # Assert first members length is odd
            if len_fst_member % 2 != 0:
                raise pytoxo.errors.ModelCSVParseError(
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
                raise pytoxo.errors.ModelCSVParseError(
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

            # Check support: only 2 different variables are supported by PyToxo
            if not 0 < len(different_variables) <= 2:
                raise pytoxo.errors.ModelCSVParseError(
                    filename,
                    "Only two variables are supported at most.",
                )
            else:
                self._variables = different_variables

            # Calculate tolerable solution error delta
            self._tolerable_solution_error_delta = (
                self.calculate_tolerable_solution_error_delta()
            )

        except IOError as e:
            raise e
        except pytoxo.errors.ModelCSVParseError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise pytoxo.errors.ModelCSVParseError(filename)

    def calculate_tolerable_solution_error_delta(self) -> float:
        """Calculates the error delta to tolerate during model resolution,
        fitted to the current model order.

        Returns
        -------
        float
            The calculated solution error delta to tolerate.
        """
        return _TOLERABLE_SOLUTION_ERROR_BASE_DELTA * pow(10, self._order)

    ########################################
    # Getters and setters for properties

    def _get_name(self) -> str:
        return self._name

    def _set_name(self, name: str) -> None:
        self._name = name

    name = property(_get_name, _set_name)

    @property
    def order(self) -> int:
        return self._order

    @property
    def penetrances(self) -> list[typing.Union[float, sympy.Expr]]:
        return self._penetrances

    @property
    def variables(self) -> list[sympy.Symbol]:
        return self._variables

    @property
    def tolerable_solution_error_delta(self) -> float:
        return self._tolerable_solution_error_delta

    ########################################

    def _max_penetrance(self) -> sympy.Expr:
        """Returns the largest of all penetrance expressions, for any real
        and positive value of the two variables and attending to the
        mathematical restrictions of the model."""

        # First, remove duplicate expressions to evaluate them
        unique_penetrances = list(set(self._penetrances))

        """The expressions must be monotonically non-decreasing and 
        sortable when `x` and `y` are real positive numbers. So, 
        to achieve the largest polynomial within the penetrance expressions,
        simply does a substitution for `x` and `y` to real positive numbers, and
        the largest numerical reduction will represent also the largest 
        symbolic expression."""
        reductions = [
            p.subs({self._variables[0]: 1, self._variables[1]: 1})
            for p in unique_penetrances
        ]  # 1 is real and positive
        return unique_penetrances[reductions.index((max(reductions)))]

    def _solve(
        self, eq_system: list[sympy.Eq], solve_timeout: typing.Union[int, bool] = True
    ) -> dict[sympy.Symbol : float]:
        """Assumes the responsibility of solving the system of equations that
        will define the values of the variables for the generation of the
        penetrance table.

        `find_max_prevalence` and `find_max_heritability` delegates here
        after conform their equations systems.

        Parameters
        ----------
        eq_system : list[sympy.Eq]
            Input equations to solve.
        solve_timeout : typing.Union[int, bool], optional (default `True`)
            A maximum timeout, as seconds, for the solver to try to resolve the
            model. If it is exceeded, the operation will be aborted. Default
            (passed as 'True') is assumed since an heuristic method related
            to the model order. Pass as 'False' or 'None' to do not use
            timeout.

        Returns
        -------
        dict[sympy.Symbol: float]
            The equation solution for the two model variables as a dict.

        Raises
        ------
        pytoxo.errors.ResolutionError
            If PyToxo is not able to solve the model.
        pytoxo.errors.UnsolvableModelError
            If the model has not solution.
        """
        # Check timeout
        if not solve_timeout:
            solve_timeout = None  # `timeout_decorator` requires exactly `None`
        elif solve_timeout == True:
            solve_timeout = (self._order + 1) ** 2 * 60  # Seconds

        # TODO: Consider add assumptions to vars real and greather than 0

        def try_to_solve(relax_dps: int = None) -> typing.Union[tuple[float], None]:
            """Tries solve the given `eq_system` equations with Sympy and
            returns solutions, filtering unreal and negative ones. The call
            to the solver is protected with a timeout to prevent
            non-termination issues of unsolvable models.

            Parameters
            ----------
            relax_dps : int (default None)
                This parameter contains the decimal digits of precision (DPS) to
                use in the internal operations delegated to MPMath library. By
                default, this library uses 15 digits. To use this value you
                could pass a 15, but better simply leave this value as a `None`
                as default. If this function does not work for you and you
                think that the model is solvable, then you could try 10 or 7,
                which have proven to be values with which the library manages
                to solve certain of the most complex models. Keep in mind that
                if you go overboard with downgrading the precision, the
                calculations could get corrupted.
            """

            @timeout_decorator.timeout_decorator.timeout(
                solve_timeout, timeout_exception=StopIteration
            )
            def solver_call():
                """Help function which encapsulate the call to the solver
                with a timeout decorator to abort futile executions of
                unsolvable models."""
                return sympy.solve(
                    eq_system,
                    self._variables[0],
                    self._variables[1],
                    manual=True,
                    rational=False,
                )

            # If applies, relax MPMath DPS to error tolerance
            if relax_dps:
                mpmath.mp.dps = relax_dps

            # Try to solve the system within the setting timeout
            try:
                sols = solver_call()
            except StopIteration:
                raise pytoxo.errors.ResolutionError(
                    cause="Exceeded timeout",
                    model_name=self.name,
                    equation=str(eq_system),
                )
            finally:
                # Re-adjust MPMath DPS
                if relax_dps:
                    mpmath.mp.dps = _MPMATH_DEFAULT_DPS

            # Discard unreal solutions
            sols = [s for s in sols if s[0].is_real and s[1].is_real]
            # Discard negative solutions
            sols = [s for s in sols if s[0] > 0 and s[1] > 0]
            # Return only one solution
            if sols:
                return sols[0]
            else:
                return None

        try:
            sol = try_to_solve()
            if not sol:
                raise pytoxo.errors.UnsolvableModelError(
                    model_name=self.name, equation=str(eq_system)
                )
        except mpmath.mp.NoConvergence:
            # Try a second time relaxing error tolerance
            try:
                sol = try_to_solve(relax_dps=10)
                if not sol:
                    raise pytoxo.errors.UnsolvableModelError(
                        model_name=self.name, equation=str(eq_system)
                    )
            except mpmath.mp.NoConvergence:
                raise pytoxo.errors.ResolutionError(
                    model_name=self.name, equation=str(eq_system)
                )

        # Errors 'UnsolvableModelError' and 'ResolutionError' are raised directly
        except pytoxo.errors.UnsolvableModelError as e:
            raise e
        except pytoxo.errors.ResolutionError as e:
            raise e

        # Wildcard drain for an unchecked situation
        except:
            raise pytoxo.errors.ResolutionError(
                model_name=self.name, equation=str(eq_system)
            )

        # Return the final achieved solution as peenetrance table object
        return {self._variables[0]: sol[0], self._variables[1]: sol[1]}

    def _check_solution(
        self, eq_system: list[sympy.Eq], sol: dict[sympy.Symbol : float]
    ) -> tuple[bool, float]:
        """Check if an achieved solution is correct.

        Parameters
        ----------
        eq_system : list[sympy.Eq]
            Equations previously solved. Right hand sides are the values
            against compare left hand sides substituted.
        sol : dict[sympy.Symbol : float]
            Achieved solution to check.

        Returns
        -------
        tuple[bool, float]
            Verification about the solution as boolean and delta desviation
            if applies.
        """
        eq1 = eq_system[0]
        eq2 = eq_system[1]

        # Substitute the solution in both equations
        eq1_lhs_eval = eq1.lhs.subs(sol).evalf()
        eq2_lhs_eval = eq2.lhs.subs(sol).evalf()

        # Calculate deltas against right hand sides
        delta1 = abs(eq1.rhs.evalf() - eq1_lhs_eval)
        delta2 = abs(eq2.rhs.evalf() - eq2_lhs_eval)

        # Check is deltas respect tolerance
        if (
            delta1 > self._tolerable_solution_error_delta
            or delta2 > self._tolerable_solution_error_delta
        ):
            return False, max(delta1, delta2)
        else:
            return True, 0

    def _build_max_prevalence_system(
        self, mafs: list[float], h: float
    ) -> list[sympy.Eq]:
        """Builds the system of equations for this model, for a maximum
        prevalence, for the given MAFs and heritability.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        h : float
            Heritability of the table.

        Returns
        -------
        list[sympy.Eq]
            System of equations to solve this model for a maximum
            prevalence, for the given MAFs and heritability.

        """
        eq1 = sympy.Eq(
            pytoxo.calculations.compute_heritability(
                self._penetrances, mafs, model_order=self._order
            ),
            h,  # It does not be casted to rational because the solver doesn't end
        )
        eq2 = sympy.Eq(self._max_penetrance(), sympy.Integer(1))
        return [eq1, eq2]  # System as a list with the equations

    def find_max_prevalence_table(
        self,
        mafs: list[float],
        h: float,
        solve_timeout: typing.Union[int, bool] = True,
        check: bool = True,
    ) -> pytoxo.ptable.PTable:
        """Computes the table whose prevalence is maximum for the given MAFs
        and heritability, and returns it within a `PTable` object.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        h : float
            Heritability of the table.
        solve_timeout : typing.Union[int, bool], optional (default True)
            A maximum timeout, as seconds, for the solver to try to resolve the
            model. If it is exceeded, the operation will be aborted. Default
            (passed as 'True') is assumed heuristically. Pass as 'False' or
            'None' to do not use timeout.
        check: bool (default True)
            Flag to control when check the solution before return it. Default
            `True`.

        Returns
        -------
        pytoxo.ptable.PTable
            Penetrance table obtained within a `PTable` object.

        Raises
        ------
        pytoxo.errors.ResolutionError
            If PyToxo is not able to solve the model.
        pytoxo.errors.UnsolvableModelError
            If the model has not solution.
        """
        eq_system = self._build_max_prevalence_system(mafs, h)

        # Delegates to the solve method
        sol = self._solve(eq_system, solve_timeout)

        # Check the solution before return the tables
        if check:
            ok, delta = self._check_solution(eq_system, sol)
            if not ok:
                raise pytoxo.errors.ResolutionError(
                    f"The calculated solution is not correct. Difference of {delta}"
                )

        # Return the final achieved solution as penetrance table object
        return pytoxo.ptable.PTable(
            model_order=self._order,
            model_penetrances=self._penetrances,
            values={
                self._variables[0]: sol[self._variables[0]],
                self._variables[1]: sol[self._variables[1]],
            },
        )

    def _build_max_heritability_system(
        self, mafs: list[float], p: float
    ) -> list[sympy.Eq]:
        """Builds the system of equations for this model, for a maximum
        heritability, for the given MAFs and prevalence.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        p : float
            Prevalence of the table.

        Returns
        -------
        list[sympy.Eq]
            System of equations to solve this model for a maximum
            heritability, for the given MAFs and prevalence.
        """
        eq1 = sympy.Eq(
            pytoxo.calculations.compute_prevalence(
                self._penetrances, mafs, model_order=self._order
            ),
            p,  # It does not be casted to rational because the solver doesn't end
        )
        eq2 = sympy.Eq(self._max_penetrance(), sympy.Integer(1))
        return [eq1, eq2]  # System as a list with the equations

    def find_max_heritability_table(
        self,
        mafs: list[float],
        p: float,
        solve_timeout: typing.Union[int, bool] = True,
        check: bool = True,
    ) -> pytoxo.ptable.PTable:
        """Computes the table whose heritability is maximum for the given MAFs
        and prevalence, and returns it within a `PTable` object.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        p : float
            Prevalence of the table.
        solve_timeout : typing.Union[int, bool], optional (default True)
            A maximum timeout, as seconds, for the solver to try to resolve the
            model. If it is exceeded, the operation will be aborted. Default
            (passed as 'True') is assumed heuristically. Pass as 'False' or
            'None' to do not use timeout.
        check: bool (default True)
            Flag to control when check the solution before return it. Default
            `True`.

        Returns
        -------
        pytoxo.ptable.PTable
            Penetrance table obtained within a `PTable` object.

        Raises
        ------
        pytoxo.errors.ResolutionError
            If PyToxo is not able to solve the model.
        pytoxo.errors.UnsolvableModelError
            If the model has not solution.
        """
        eq_system = self._build_max_heritability_system(mafs, p)

        # Delegates to the solve method
        sol = self._solve(eq_system, solve_timeout)

        # Check the solution before return the tables
        if check:
            ok, delta = self._check_solution(eq_system, sol)
            if not ok:
                raise pytoxo.errors.ResolutionError(
                    f"The calculated solution is not correct. Difference of {delta}"
                )

        # Return the final achieved solution as penetrance table object
        return pytoxo.ptable.PTable(
            model_order=self._order,
            model_penetrances=self._penetrances,
            values={
                self._variables[0]: sol[self._variables[0]],
                self._variables[1]: sol[self._variables[1]],
            },
        )
