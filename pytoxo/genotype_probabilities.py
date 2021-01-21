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

"""Genotype probabilities operation."""

import functools
import itertools
import operator

from sympy import sympify


def genotype_probabilities(mafs: [float]) -> [float]:
    """Compute the probabilities associated with all genotype combinations given
    each MAF (minor allele frequency).

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
    mafs : [float]
        Minor allele frequencies array.

    Returns
    -------
    [float]
        Array with the probabilities of all possible allele combinations.
    """
    af_zip = []  # Zip with pairs of alleles frequencies as `(m, M)` tuples
    for maf in mafs:
        m = sympify(maf)
        M = 1 - sympify(maf)
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
    for o in gen_probs_cp:
        gen_probs_cp_red.append(functools.reduce(operator.mul, o, 1))

    return gen_probs_cp_red
