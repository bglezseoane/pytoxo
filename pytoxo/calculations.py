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

"""PyToxo calculations module: some math-based methods used in different
parts of the library.

If possible, internal Sympy classes and methods are chosen instead of Python
built-ins for optimization.
"""

import functools
import inspect
import itertools
from typing import List, Union

from sympy import Integer, Rational, Expr, Add, Mul, Pow, nsimplify, simplify

import pytoxo.util
from pytoxo.errors import GenericCalculationError


def genotype_probabilities(
    mafs: Union[List[Rational], List[float]], model_order: int = None
) -> Union[List[Rational], List[Expr]]:
    """Computes the probabilities associated with all genotype combinations
    given each MAF (minor allele frequency).

    The genotype frequency can be calculated from the allele frequencies. It
    is usually assumed for its calculation that the Hardy-Weinberg equilibrium
    is fulfilled, that is, if the frequency of allele `A` is `M` and the
    frequency of allele `a` is `m`, the probability of having the `AA` genotype
    is `M^2`, having `Aa` is `2 * M * m` and having `aa` is `m^2` (this could be
    extended to any order of interaction, including the probabilities of each
    allele in the products).

    Therefore, given the MAF (minor allele frequency), the frequencies of all
    possible genotypes can be calculated. For example, if MAF is 0.25, then
    `M = 0.25`, and `m = 1 - 0.25 = 0.75`. From there, the probabilities of
    all possible allele combinations (of all possible genotypes) can be
    calculated.

    Parameters
    ----------
    mafs : Union[list[Rational], list[float]]
        Minor allele frequencies array. Accept rationals and floats.
    model_order : int, optional
        The order of the implicated model. This optional parameter is used to
        optimize the computation. Default is `None` and supposes a generalist
        optimization strategy, less optimal.

    Returns
    -------
    Union[list[Rational], list[Expr]]
        Array with the probabilities of all possible allele combinations as
        rational numbers.

    Raises
    ------
    GenericCalculationError
        On any error situation.
    """
    try:
        mafs = [nsimplify(m) for m in mafs]  # Assert rationals

        af_zip = []  # Zip with pairs of alleles frequencies as `(m, M)` tuples
        for maf in mafs:
            m = maf
            M = Add(Integer(1), Mul(Integer(-1), maf))
            af_zip.append((m, M))

        gen_probs = (
            []
        )  # Genotype probabilities as `[M ** 2, 2 * M * m, m ** 2]` sub-lists
        for allele_pair in af_zip:
            m = allele_pair[0]
            M = allele_pair[1]
            gen_probs.append(
                [Pow(M, Integer(2)), Mul(Integer(2), M, m), Pow(m, Integer(2))]
            )

        # Cartesian product of the probabilities
        gen_probs_cps = list(itertools.product(*gen_probs))

        # Reduction of the cartesian product with a multiplication for each sub-list
        gen_probs_cp_red = []
        for gen_probs_cp in gen_probs_cps:
            gen_probs_cp_red.append(functools.reduce(Mul, gen_probs_cp, Integer(1)))

        return gen_probs_cp_red

    except:
        frame = inspect.currentframe()
        raise GenericCalculationError(inspect.getframeinfo(frame).function)


def compute_prevalence(
    penetrances: Union[List[Expr], List[Rational], List[float]],
    mafs: Union[List[Rational], List[float]] = None,
    gp: Union[List[Rational], List[float]] = None,
    model_order: int = None,
) -> Union[Rational, Expr]:
    """Tries to compute the prevalence for a given penetrance list.

    One of `mafs` and `gp` is required. Is the both are used, `mafs` is ignored.

    Parameters
    ----------
    penetrances : Union[list[Expr], list[Rational], list[float]]
        Penetrance values array.
    mafs : Union[list[Rational], list[float]], optional
        Minor allele frequencies array. Ignored if `gp` is used.
    gp : Union[list[Rational], list[float]], optional
        Genotype probabilities array.
    model_order : int, optional
        The order of the implicated model. This optional parameter is used to
        optimize the computation. Default is `None` and supposes a generalist
        optimization strategy, less optimal.

    Returns
    -------
    Union[Rational, Expr]
        Prevalence of the penetrance table. Returns a rational if it is
        possible to solve the expression numerically and the expression if not.

    Raises
    ------
    ValueError
        If none of `mafs` and `gp` are passed as parameter.

    Raises
    ------
    GenericCalculationError
        On any error situation.
    """
    try:
        if not gp:
            if not mafs:
                raise ValueError(
                    "One of `mafs` and `gp` is required. Is the both are used, `mafs` is ignored."
                )
            else:
                gp = genotype_probabilities(mafs)
        else:
            gp = [nsimplify(p) for p in gp]  # Assert rationals

        penetrances = [nsimplify(p) for p in penetrances]  # Assert rationals

        # `penetrances .* gp`
        prods = []
        for pen, prob in zip(penetrances, gp):
            prods.append(Mul(pen, prob))

        addition = Add(*prods)  # Addition of all the elements of the product array

        """Try to solve. Innocuous if the expression is still symbolic and saves a
        lot of time if we are dealing with an expression that is already reducible 
        to a numeric value"""
        addition = addition.evalf().nsimplify()

        return _try_to_simplify(
            addition, model_order
        )  # Return simplified to optimize next steps

    except:
        frame = inspect.currentframe()
        raise GenericCalculationError(inspect.getframeinfo(frame).function)


def compute_heritability(
    penetrances: Union[List[Expr], List[Rational], List[float]],
    mafs: Union[List[Rational], List[float]],
    model_order: int = None,
) -> Union[Rational, Expr]:
    """Tries to compute the heritability for a given penetrance table
    defined by its values.

    Parameters
    ----------
    penetrances : Union[list[Expr], list[Rational], list[float]]
        Penetrance values array.
    mafs : Union[list[Rational], list[float]], optional
        Minor allele frequencies array.
    model_order : int, optional
        The order of the implicated model. This optional parameter is used to
        optimize the computation. Default is `None` and supposes a generalist
        optimization strategy, less optimal.

    Returns
    -------
    Union[Rational, Expr]
        Heritability of the penetrance table. Returns a rational if it is
        possible to solve the expression numerically and the expression if not.

    Raises
    ------
    GenericCalculationError
        On any error situation.
    """
    try:
        penetrances = [nsimplify(p) for p in penetrances]  # Assert rationals

        gp = genotype_probabilities(mafs)
        p = compute_prevalence(penetrances, mafs, gp)

        # `(penetrances - p).^2 .* gp`
        prods = []
        for pen, prob in zip(penetrances, gp):
            prods.append(Mul(Pow(Add(pen, Mul(Integer(-1), p)), Integer(2)), prob))

        # Denominator of the final expression
        denom = _try_to_simplify(
            Pow(Mul(p, Add(Integer(1), Mul(Integer(-1), p))), Integer(-1), model_order)
        )

        # `prods / denom`, because `denom` is a negative pow
        mult = Mul(Add(*prods), denom)

        """Try to solve. Innocuous if the expression is still symbolic and saves a
        lot of time if we are dealing with an expression that is already reducible 
        to a numeric value"""
        mult = mult.evalf().nsimplify()

        return _try_to_simplify(
            mult, model_order
        )  # Return simplified to optimize next steps

    except:
        frame = inspect.currentframe()
        raise GenericCalculationError(inspect.getframeinfo(frame).function)


def _try_to_simplify(expr: Expr, timeout: int = 60, model_order: int = None):
    """Help function which encapsulate the a call to an attempt to simplify an
    expression with a maximum timeout. If the expression could not be simplified
    in that time, it returns it as it was at the beginning, else return it
    simplified. In this way it is ensured that the simplifications of this
    module help the efficiency of the program and not the opposite. It also
    prevents the program from seemingly stuck trying to simplify overly
    complex expressions.

    Parameters
    ----------
    expr : Expr
        Expression to simplify.
    timeout : int, optional (default 60)
        Maximum timeout for the try, as seconds. Default 1 minute as 60 seconds.
        If the `model_order` parameter is used, a better value for this is
        also achieved. Use `None` to don't use a timeout.
    model_order : int, optional
        The order of the implicated model. This optional parameter is used to
        optimize the computation. Default is `None` and supposes a generalist
        optimization strategy, less optimal.

    Returns
    -------
    Expr
        The input expression, simplified or not.
    """
    if model_order:
        timeout = model_order * 60  # Seconds

    @pytoxo.util.timeout(timeout)
    def simplify_call():
        return simplify(expr)

    try:
        return simplify_call()
    except TimeoutError:
        return expr
