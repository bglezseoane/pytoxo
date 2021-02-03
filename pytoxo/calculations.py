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
    penetrances: list[float], mafs: list[float], gp: list[float] = None
) -> typing.Union[float, sympy.Expr]:
    """Tries to compute the prevalence for a given penetrance table
    defined by its values.

    Parameters
    ----------
    penetrances : list[float]
        Penetrance values array.
    mafs : list[float]
        Minor allele frequencies array. Ignored if `gp` is used.
    gp : list[float], optional
        Genotype probabilities array.

    Returns
    -------
    typing.Union[float, sympy.Expr]
        Prevalence of the penetrance table. Returns a float if the expression is
        already solved and the expression itself as a Sympy expression
        object if not.
    """
    if not gp:
        gp = genotype_probabilities(mafs)

    res = np.sum(np.multiply(penetrances, gp))

    if type(res) is np.ndarray:
        res = float(res)

    return res


def compute_heritability(
    penetrances: list[float], mafs: list[float]
) -> typing.Union[float, sympy.Expr]:
    """Tries to compute the heritability for a given penetrance table
    defined by its values.

    Parameters
    ----------
    penetrances : list[float]
        Penetrance values array.
    mafs : list[float]
        Minor allele frequencies array.

    Returns
    -------
    typing.Union[float, sympy.Expr]
        Heritability of the penetrance table. Returns a float if the expression
        is already solved and the expression itself as a Sympy expression
        object if not.
    """
    gp = genotype_probabilities(mafs)
    p = compute_prevalence(penetrances, mafs, gp)

    res = np.sum(np.multiply(np.power(np.subtract(penetrances, p), 2), gp)) / (
        p * (1 - p)
    )

    if type(res) is np.ndarray:
        res = float(res)

    return res
