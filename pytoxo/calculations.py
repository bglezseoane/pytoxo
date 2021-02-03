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
parts of the library."""

import functools
import itertools
import operator
import typing

import sympy
import numpy as np


def genotype_probabilities(
    mafs: typing.Union[list[sympy.Rational], list[float]]
) -> list[sympy.Rational]:
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
    mafs : typing.Union[list[sympy.Rational], list[float]]
        Minor allele frequencies array. Accept rationals and floats.

    Returns
    -------
    list[sympy.Rational]
        Array with the probabilities of all possible allele combinations as
        rational numbers.
    """
    mafs = [sympy.nsimplify(m) for m in mafs]  # Assert rationals

    af_zip = []  # Zip with pairs of alleles frequencies as `(m, M)` tuples
    for maf in mafs:
        m = sympy.sympify(maf)
        M = 1 - sympy.sympify(maf)
        af_zip.append((m, M))

    gen_probs = []  # Genotype probabilities as `[M ** 2, 2 * M * m, m ** 2]` sub-lists
    for allele_pair in af_zip:
        m = allele_pair[0]
        M = allele_pair[1]
        gen_probs.append([M ** 2, 2 * M * m, m ** 2])

    # Cartesian product of the probabilities
    gen_probs_cp = list(itertools.product(*gen_probs))

    # Reduction of the cartesian product with a multiplication for each sub-list
    gen_probs_cp_red = []
    for sl in gen_probs_cp:
        gen_probs_cp_red.append(functools.reduce(operator.mul, sl, 1))

    return gen_probs_cp_red


def compute_prevalence(
    penetrances: typing.Union[list[sympy.Expr], list[sympy.Rational], list[float]],
    mafs: typing.Union[list[sympy.Rational], list[float]],
    gp: typing.Union[list[sympy.Rational], list[float]] = None,
) -> typing.Union[sympy.Rational, sympy.Expr]:
    """Tries to compute the prevalence for a given penetrance list.

    One of `mafs` and `gp` is required. Is the both are used, `mafs` is ignored.

    Parameters
    ----------
    penetrances : typing.Union[list[sympy.Expr], list[sympy.Rational], list[float]]
        Penetrance values array.
    mafs : typing.Union[list[sympy.Rational], list[float]], optional
        Minor allele frequencies array. Ignored if `gp` is used.
    gp : typing.Union[list[sympy.Rational], list[float]], optional
        Genotype probabilities array.

    Returns
    -------
    typing.Union[sympy.Rational, sympy.Expr]
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
        gp = [sympy.nsimplify(p) for p in gp]  # Assert rationals

    penetrances = [sympy.nsimplify(p) for p in penetrances]  # Assert rationals

    # Product of each value of the two arrays
    prods = []
    for pen, prob in zip(penetrances, gp):
        prods.append(pen * prob)

    return sympy.Add(
        *prods
    )  # Return the addition of all the elements of the product array


def compute_heritability(
    penetrances: typing.Union[list[sympy.Expr], list[sympy.Rational], list[float]],
    mafs: typing.Union[list[sympy.Rational], list[float]],
) -> typing.Union[sympy.Rational, sympy.Expr]:
    """Tries to compute the heritability for a given penetrance table
    defined by its values.

    Parameters
    ----------
    penetrances : typing.Union[list[sympy.Expr], list[sympy.Rational], list[float]]
        Penetrance values array.
    mafs : typing.Union[list[sympy.Rational], list[float]], optional
        Minor allele frequencies array.

    Returns
    -------
    typing.Union[sympy.Rational, sympy.Expr]
        Heritability of the penetrance table. Returns a rational if it is
        possible to solve the expression numerically and the expression if not.
    """
    penetrances = [sympy.nsimplify(p) for p in penetrances]  # Assert rationals

    gp = genotype_probabilities(mafs)
    p = compute_prevalence(penetrances, mafs, gp)

    sub_squares = []
    # Subtraction squares of each value of the two arrays
    for pen in penetrances:
        sub_squares.append(sympy.Pow(sympy.Add(pen, sympy.Mul(-1, p)), 2))

    # Product of each value of the two arrays
    prods = []
    for ssq, prob in zip(sub_squares, gp):
        prods.append(ssq * prob)

    # Denominator of the final expression
    denom = sympy.Pow(sympy.Mul(p, sympy.Add(1, sympy.Mul(p, -1))), -1)

    return sympy.Mul(sympy.Add(*prods), denom)
