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

"""Epistasis model definition."""

import itertools
import os
import string
from typing import Dict, List, Tuple, Union

import mpmath
import numpy
import sympy

import pytoxo.calculations
import pytoxo.errors
import pytoxo.ptable
import pytoxo.util

_MPMATH_DEFAULT_DPS = 15
_TOLERABLE_SOLUTION_ERROR_BASE_DELTA = 1e-16  # It is fitted then to model's order
_TOLERABLE_SOLUTION_ERROR_MAX_DELTA = 1e-8


class Model:
    """Representation of an epistasis model."""

    def __init__(
        self,
        filename: str = None,
        genotypes_dict: Dict[str, str] = None,
        definitions: Union[List[str], numpy.array] = None,
        probabilities: Union[List[str], numpy.array] = None,
        model_name: str = None,
    ):
        """Initializes a model object.

        This constructor can be used in three modes:

        1. Generating the model object from its representation as a CSV file.
        2. Generating the model object from a dictionary with the genotype
            definitions and the expressions of the probabilities of the
            genotypes (`genotypes_dict`).
        3. Generating the model object from two separated lists or Numpy
            arrays with the genotype definitions and the expressions of the
            probabilities of the genotypes (`definitions` and `probabilities`).

        Generating from its representation as a CSV file
        ------------------------------------------------

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
        file. Math syntax should be the used by Sympy.

        Empty lines, as well as lines starting with '#' will be ignored.

        To use this constructor method, use the `filename` parameter
        and let others as default `None` but optional `model_name`.


        Generating from a dictionary of genotype definitions and probabilities
        ----------------------------------------------------------------------

        Uses a Python dictionary, `genotypes_dict`, with the genotypes
        definitions and its associated probabilities, both as strings. E.g.
        of a member of the dict: `"AABBCc": "x*(1+y)"`.

        Probability is expressed as a function of two variables, represented
        as alphabetical characters. Only two variables are supported at most.
        Any pair of alphabetic characters can be used to represent these
        variables, but they must be the same two names for the entire model
        file. Math syntax should be the used by Sympy.

        To use this constructor method, use the `genotypes_dict` parameter
        and let others as default `None` but optional `model_name`.


        Generating from two separated Python lists or Numpy arrays
        ----------------------------------------------------------

        Uses two separated Python lists or Numpy arrays, `definitions` and
        `probabilities`, with the genotypes definitions and its associated
        probabilities, respectively, both as strings. The link between the
        members of both sets is determined by the order. Examples of the two
        sets, respectively: "AABBCc" and "x*(1+y)".

        Probability is expressed as a function of two variables, represented
        as alphabetical characters. Only two variables are supported at most.
        Any pair of alphabetic characters can be used to represent these
        variables, but they must be the same two names for the entire model
        file. Math syntax should be the used by Sympy.

        To use this constructor method, use the `definitions` and
        `probabilities` parameters and let others as default `None` but
        optional `model_name`.


        Parameters
        ----------
        filename : str, optional
            The path of the text file with the model.
        genotypes_dict : dict[str, str], optional
            Dict with genotypes definitions and its associated probabilities.
        definitions: Union[list[str], numpy.array], optional
            List ot Numpy array with the genotype definitions as strings.
            Should be used with `probabilities` parameter.
        probabilities: Union[list[str], numpy.array], optional
            List ot Numpy array with the genotype probabilities as strings.
            Should be used with `definitions` parameter.
        model_name: str, optional
            The name to identify the model. Optional. Using the `filename`
            can be automatically deduced attending to the file name, if this
            parameter is not used. Using the `genotypes_dict` constructor and
            letting this parameter unused the result name would be "unnamed".

        Raises
        ------
        pytoxo.errors.BadFormedModelError
            If the parse tentative fails due to the file is not well formed.
        ValueError
            Bad parameters combination.
        IOError
           If the input file is not found or the reading tentative fails due
           to other unexpected operative system level cause.
        """
        # Parse parameters use to determine constructor way to use. Ignore
        # for now `model_name`
        exclusive_set_params = [
            True
            for p in [filename, genotypes_dict, definitions, probabilities]
            if p is not None
        ]
        if len(exclusive_set_params) == 1 and filename:
            constructor_mode = 1
        elif len(exclusive_set_params) == 1 and genotypes_dict:
            constructor_mode = 2
        elif (
            len(exclusive_set_params) == 2
            and definitions
            is not None  # Need to be explicit because could be Numpy array
            and probabilities
            is not None  # Need to be explicit because could be Numpy array
        ):
            constructor_mode = 3
        else:
            raise ValueError("Bad parameters use. Revise documentation.")

        # Basic skeleton
        self._name = None
        self._order = None  # Number of loci defined in the model
        self._penetrances = (
            []
        )  # List of symbolic expressions representing the epistatic model
        self._variables = (
            []
        )  # List of symbolic variable names used throughout the model

        # Check custom model name
        if model_name:
            self._name = str(model_name)

        # Delegates object composition attending to the used constructor way
        if constructor_mode == 1:
            self._parse_model_file(filename)
        elif constructor_mode == 2:
            self._parse_genotypes_dict(genotypes_dict)
        elif constructor_mode == 3:
            self._parse_genotypes_sets(definitions, probabilities)

    def _parse_model_file(self, filename: str) -> None:
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
        file. Math syntax should be the used by Sympy.

        Empty lines, as well as lines starting with '#' will be ignored.

        Parameters
        ----------
        filename : str
            The path of the text file with the model.

        Raises
        ------
        pytoxo.errors.BadFormedModelError
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
                raise pytoxo.errors.BadFormedModelError(
                    filename, "File without content."
                )

            # Split lines around the comma (',')
            fst_members = [line.split(",")[0].strip() for line in lines]
            snd_members = [line.split(",")[1].strip() for line in lines]

            # Save the name of the model, if a custom one is not used
            if not self._name:
                self._name = os.path.basename(filename).split(".")[0]

            # Delegate to the helper
            self._parse_helper(
                genotypes=fst_members,
                probabilities=snd_members,
                exception_object_to_raise=filename,
            )

        except IOError as e:
            raise e
        except pytoxo.errors.BadFormedModelError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise pytoxo.errors.BadFormedModelError(filename)

    def _parse_genotypes_dict(self, genotypes_dict: Dict[str, str]) -> None:
        """Takes the responsibility of the initializer to parse the
        genotypes dict.

        Uses a Python dictionary, `genotypes_dict`, with the genotypes
        definitions and its associated probabilities, both as strings. E.g.
        of a member of the dict: `"AABBCc": "x*(1+y)"`.

        Probability is expressed as a function of two variables, represented
        as alphabetical characters. Only two variables are supported at most.
        Any pair of alphabetic characters can be used to represent these
        variables, but they must be the same two names for the entire model
        file. Math syntax should be the used by Sympy.

        Parameters
        ----------
        genotypes_dict : dict[str, str]
            Dict with genotypes definitions and its associated probabilities.

        Raises
        ------
        pytoxo.errors.BadFormedModelError
            If the input genotypes dict is bad-formed.
        """
        try:
            # Delegate to the helper
            self._parse_helper(
                genotypes=list(genotypes_dict.keys()),
                probabilities=list(genotypes_dict.values()),
                exception_object_to_raise=genotypes_dict,
            )

            # Revise unnamed models
            if not self._name:
                self._name = "unnamed"

        except pytoxo.errors.BadFormedModelError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise pytoxo.errors.BadFormedModelError(genotypes_dict)

    def _parse_genotypes_sets(
        self,
        definitions: Union[List[str], numpy.array],
        probabilities: Union[List[str], numpy.array],
    ) -> None:
        """Takes the responsibility of the initializer to parse the
        two sets with genotype definitions and probabilities.

        Uses two separated Python lists or Numpy arrays, `definitions` and
        `probabilities`, with the genotypes definitions and its associated
        probabilities, respectively, both as strings. The link between the
        members of both sets is determined by the order. Examples of the two
        sets, respectively: "AABBCc" and "x*(1+y)".

        Probability is expressed as a function of two variables, represented
        as alphabetical characters. Only two variables are supported at most.
        Any pair of alphabetic characters can be used to represent these
        variables, but they must be the same two names for the entire model
        file. Math syntax should be the used by Sympy.

        To use this constructor method, use the `definitions` and
        `probabilities` parameters and let others as default `None` but
        optional `model_name`.


        Parameters
        ----------
        definitions: Union[list[str], numpy.array]
            List ot Numpy array with the genotype definitions as strings.
            Should be used with `probabilities` parameter.
        probabilities: Union[list[str], numpy.array]
            List ot Numpy array with the genotype probabilities as strings.
            Should be used with `definitions` parameter.

        Raises
        ------
        pytoxo.errors.BadFormedModelError
            If the input genotype sets are bad-formed.
        """
        try:
            # Transform Numpy arrays to lists
            if type(definitions) == numpy.array:
                definitions = definitions.tolist()
            if type(probabilities) == numpy.array:
                probabilities = probabilities.tolist()

            # Delegate to the helper
            self._parse_helper(
                genotypes=list(definitions),
                probabilities=list(probabilities),
                exception_object_to_raise=(definitions, probabilities),
            )

            # Revise unnamed models
            if not self._name:
                self._name = "unnamed"

        except pytoxo.errors.BadFormedModelError as e:
            raise e
        except:  # Generic drain for unchecked parsing errors
            raise pytoxo.errors.BadFormedModelError((definitions, probabilities))

    def _parse_helper(
        self,
        genotypes: List[str],
        probabilities: List[str],
        exception_object_to_raise: Union[str, Dict[str, str], Tuple[list, list]],
    ) -> None:
        """Complete the parse process doing the needed checks and
        transformations.

        Parameters
        ----------
        genotypes: list[str]
            List with the genotype definitions as strings.
        probabilities: list[str]
            List with the genotype associated probability expressions as strings.
        exception_object_to_raise: Union[str, dict[str, str], tuple[list, list]]
            Object to add to the exception body.

        Raises
        ------
        pytoxo.errors.BadFormedModelError
            On any fail.
        """
        # Catch the order of the model
        len_genotype = len(genotypes[0])
        # Assert all the first members have the same length
        for genotype in genotypes:
            if len(genotype) != len_genotype:
                raise pytoxo.errors.BadFormedModelError(
                    exception_object_to_raise,
                    "All the genotypes should have the same order.",
                )
        # Assert first members length is odd
        if len_genotype % 2 != 0:
            raise pytoxo.errors.BadFormedModelError(
                exception_object_to_raise, "Bad genotype specification."
            )
        # Save the order
        self._order = len_genotype // 2

        """Sort probabilities attending to alphabetical sort of associated 
        genotypes. Capital letters first. This is necessary to assert the 
        association between genotype definitions and probabilities during the 
        calculus process and in the final penetrance table"""
        genotypes_probabilities = [(g, p) for g, p in zip(genotypes, probabilities)]
        genotypes_probabilities.sort(key=lambda i: i[0])  # Sort attending to genotypes
        probabilities_sorted = [p for (_, p) in genotypes_probabilities]
        # Save the penetrances as symbolic expressions
        try:
            self._penetrances = [
                sympy.sympify(probability) for probability in probabilities_sorted
            ]
        except sympy.SympifyError:
            raise pytoxo.errors.BadFormedModelError(
                exception_object_to_raise, "Bad probability expression syntax."
            )

        # Save the variables of the used expressions
        all_variables = []
        for probability in probabilities_sorted:
            all_variables.append(sympy.symbols([i for i in probability if i.isalpha()]))

        # Reduce the variables to avoid replication
        all_variables = [item for subl in all_variables for item in subl]  # Flat list
        different_variables = []
        for variable in all_variables:
            """Not use other more direct method like
            `list(set(all_variables))` because the variable apparition
            order must be preserved"""
            if variable not in different_variables:
                different_variables.append(variable)

        # Check support: only 2 different variables are supported by PyToxo
        if not 0 < len(different_variables) <= 2:
            raise pytoxo.errors.BadFormedModelError(
                exception_object_to_raise,
                "Only two variables are supported at most.",
            )
        else:
            self._variables = different_variables

        # Calculate tolerable solution error delta
        self._tolerable_solution_error_delta = (
            self.calculate_tolerable_solution_error_delta()
        )

    def calculate_tolerable_solution_error_delta(self) -> float:
        """Calculates the error delta to be tolerated during the model
        resolution, adjusted to the current model order. If the adjusted
        delta is greater than the maximum delta, returns the maximum delta.
        The latter avoid absurd delta error for very large models.

        Uses the constants `_TOLERABLE_SOLUTION_ERROR_BASE_DELTA` and
        `_TOLERABLE_SOLUTION_ERROR_MAX_DELTA` to calculate the fitted delta.

        Returns
        -------
        float
            The resolution error delta to tolerate.
        """
        return min(
            _TOLERABLE_SOLUTION_ERROR_BASE_DELTA * pow(10, self._order),
            _TOLERABLE_SOLUTION_ERROR_MAX_DELTA,
        )

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
    def penetrances(self) -> List[Union[float, sympy.Expr]]:
        return self._penetrances

    @property
    def variables(self) -> List[sympy.Symbol]:
        return self._variables

    @property
    def tolerable_solution_error_delta(self) -> float:
        return self._tolerable_solution_error_delta

    ########################################

    def __hash__(self):
        return hash(
            hash(self._name)
            + hash(self._order)
            + hash(str(self._penetrances))
            + hash(str(self._variables))
            + hash(str(self._tolerable_solution_error_delta))
        )

    def __eq__(self, other):
        return hash(self) == hash(other)

    def calculate_genotypes(self) -> List[str]:
        """Calculate the list of genotypes to the given model attending to
        the model order.

        Uses default alphabetical sort. Capital letters first.

        Returns
        -------
        list[str]
            List of genotypes.
        """
        # Deduce alleles attending to the table order
        # TODO: Is `len(ascii_lowercase) < self._order` possible?
        letters = list(string.ascii_lowercase[: self._order])
        alleles = []
        for letter in letters:
            lower = letter
            upper = letter.upper()
            alleles.append([f"{upper}{upper}", f"{upper}{lower}", f"{lower}{lower}"])
        # Generate genotypes tracing the alleles cartesian product
        return [str("".join(p)) for p in list(itertools.product(*alleles))]

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
        self, eq_system: List[sympy.Eq], solve_timeout: Union[int, bool] = True
    ) -> Dict[sympy.Symbol, float]:
        """Assumes the responsibility of solving the system of equations that
        will define the values of the variables for the generation of the
        penetrance table.

        `find_max_prevalence` and `find_max_heritability` delegates here
        after conform their equations systems.

        Parameters
        ----------
        eq_system : list[sympy.Eq]
            Input equations to solve.
        solve_timeout : Union[int, bool], optional (default `True`)
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
            solve_timeout = None  # Required exactly `None`
        elif type(solve_timeout) == bool and solve_timeout == True:
            solve_timeout = (self._order + 1) ** 2 * 60  # Seconds

        # TODO: Consider add assumptions to vars real and greather than 0

        def try_to_solve(relax_dps: int = None) -> Union[Tuple[float], None]:
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

            @pytoxo.util.timeout(solve_timeout)
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
            except TimeoutError:
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
        self, eq_system: List[sympy.Eq], sol: Dict[sympy.Symbol, float]
    ) -> Tuple[bool, float]:
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
        self, mafs: List[float], h: float
    ) -> List[sympy.Eq]:
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
        mafs: List[float],
        h: float,
        solve_timeout: Union[int, bool] = True,
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
        solve_timeout : Union[int, bool], optional (default True)
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
        ValueError
            With a bad parameter configuration.
        """
        self.check_find_table_parameters(mafs, h, solve_timeout, check)

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
            model_genotypes=self.calculate_genotypes(),
            model_penetrances=self._penetrances,
            values={
                self._variables[0]: sol[self._variables[0]],
                self._variables[1]: sol[self._variables[1]],
            },
            mafs=mafs,
        )

    def _build_max_heritability_system(
        self, mafs: List[float], p: float
    ) -> List[sympy.Eq]:
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
        mafs: List[float],
        p: float,
        solve_timeout: Union[int, bool] = True,
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
        solve_timeout : Union[int, bool], optional (default True)
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
        ValueError
            With a bad parameter configuration.
        """
        self.check_find_table_parameters(mafs, p, solve_timeout, check)

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
            model_genotypes=self.calculate_genotypes(),
            model_penetrances=self._penetrances,
            values={
                self._variables[0]: sol[self._variables[0]],
                self._variables[1]: sol[self._variables[1]],
            },
            mafs=mafs,
        )

    def check_find_table_parameters(
        self,
        mafs: List[float],
        h_or_p: float,
        solve_timeout: Union[int, bool] = True,
        check: bool = True,
    ) -> None:
        """Check the input parameters for a find table operation.

        Parameters
        ----------
        mafs : list[float]
            Array of floats representing the MAF of each locus.
        h_or_p : float
            Heritability or prevalence of the table.
        solve_timeout : Union[int, bool], optional (default True)
            A maximum timeout, as seconds, for the solver to try to resolve the
            model. If it is exceeded, the operation will be aborted. Default
            (passed as 'True') is assumed heuristically. Pass as 'False' or
            'None' to do not use timeout.
        check: bool (default True)
            Flag to control when check the solution before return it. Default
            `True`.

        Raises
        ------
        ValueError
            With a bad parameter configuration.
        """
        if type(mafs) is not list:
            raise ValueError("MAFs should be a list of floats.")
        for maf in mafs:
            if type(maf) is not float:
                raise ValueError("MAFs should be a list of floats.")
        if len(mafs) != self._order:
            raise ValueError("The MAFs have to be as many as the order has the model.")
        for maf in mafs:
            if not 0 <= maf <= 0.5:
                raise ValueError("MAFs should be a probability between 0 and 0.5.")

        if type(h_or_p) is not float:
            raise ValueError("Heritability or prevalence should be a float.")
        if not 0 <= h_or_p <= 1:
            raise ValueError(
                "Heritability or prevalence should be probability. I.e. a float between 0 and 1."
            )

        if type(solve_timeout) is not int and type(solve_timeout) is not bool:
            raise ValueError("Timeout parameter should be an integer or a boolean.")
        if type(solve_timeout) is int and solve_timeout < 0:
            raise ValueError("Timeout should be a positive integer.")

        if type(check) is not bool:
            raise ValueError("The check flag should be a boolean.")
