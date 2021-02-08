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

"""PyToxo calculations module: some math-based methods used in different
parts of the library.

If possible, internal Sympy classes and methods are chosen instead of Python
built-ins for optimization.
"""

import functools
import itertools
import typing

import timeout_decorator
from sympy import Integer, Rational, Expr, Add, Mul, Pow, nsimplify, simplify


def genotype_probabilities(
    mafs: typing.Union[list[Rational], list[float]]
) -> typing.Union[list[Rational], list[Expr]]:
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
    mafs : typing.Union[list[Rational], list[float]]
        Minor allele frequencies array. Accept rationals and floats.

    Returns
    -------
    typing.Union[list[Rational], list[Expr]]
        Array with the probabilities of all possible allele combinations as
        rational numbers.
    """
    mafs = [nsimplify(m) for m in mafs]  # Assert rationals

    af_zip = []  # Zip with pairs of alleles frequencies as `(m, M)` tuples
    for maf in mafs:
        m = maf
        M = Add(Integer(1), Mul(Integer(-1), maf))
        af_zip.append((m, M))

    gen_probs = []  # Genotype probabilities as `[M ** 2, 2 * M * m, m ** 2]` sub-lists
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


def compute_prevalence(
    penetrances: typing.Union[list[Expr], list[Rational], list[float]],
    mafs: typing.Union[list[Rational], list[float]] = None,
    gp: typing.Union[list[Rational], list[float]] = None,
) -> typing.Union[Rational, Expr]:
    """Tries to compute the prevalence for a given penetrance list.

    One of `mafs` and `gp` is required. Is the both are used, `mafs` is ignored.

    Parameters
    ----------
    penetrances : typing.Union[list[Expr], list[Rational], list[float]]
        Penetrance values array.
    mafs : typing.Union[list[Rational], list[float]], optional
        Minor allele frequencies array. Ignored if `gp` is used.
    gp : typing.Union[list[Rational], list[float]], optional
        Genotype probabilities array.

    Returns
    -------
    typing.Union[Rational, Expr]
        Prevalence of the penetrance table. Returns a rational if it is
        possible to solve the expression numerically and the expression if not.

    Raises
    ------
    ValueError
        If none of `mafs` and `gp` are passed as parameter.
    """
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

    return _try_to_simplify(addition)  # Return simplified to optimize next steps


def compute_heritability(
    penetrances: typing.Union[list[Expr], list[Rational], list[float]],
    mafs: typing.Union[list[Rational], list[float]],
) -> typing.Union[Rational, Expr]:
    """Tries to compute the heritability for a given penetrance table
    defined by its values.

    Parameters
    ----------
    penetrances : typing.Union[list[Expr], list[Rational], list[float]]
        Penetrance values array.
    mafs : typing.Union[list[Rational], list[float]], optional
        Minor allele frequencies array.

    Returns
    -------
    typing.Union[Rational, Expr]
        Heritability of the penetrance table. Returns a rational if it is
        possible to solve the expression numerically and the expression if not.
    """
    penetrances = [nsimplify(p) for p in penetrances]  # Assert rationals

    gp = genotype_probabilities(mafs)
    p = compute_prevalence(penetrances, mafs, gp)

    # `(penetrances - p).^2 .* gp`
    prods = []
    for pen, prob in zip(penetrances, gp):
        prods.append(
            _try_to_simplify(Mul(Pow(Add(pen, Mul(Integer(-1), p)), Integer(2)), prob))
        )

    # Denominator of the final expression
    denom = _try_to_simplify(
        Pow(Mul(p, Add(Integer(1), Mul(Integer(-1), p))), Integer(-1))
    )

    # `prods / denom`, because `denom` is a negative pow
    mult = Mul(Add(*prods), denom)

    a = _try_to_simplify(mult)  # Return simplified to optimize next steps
    return a


def _try_to_simplify(expr: Expr, timeout: int = 30):
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

    Returns
    -------
    Expr
        The input expression, simplified or not.

    """

    @timeout_decorator.timeout_decorator.timeout(
        timeout, timeout_exception=StopIteration
    )
    def simplify_call():
        return simplify(expr)

    try:
        return simplify_call()
    except StopIteration:
        return expr
