"""Balancing chemical equations and reaction stoichiometry.

Balancing finds the integer null space of the element-composition matrix
using exact rational Gauss-Jordan elimination, so there's no floating-point
round-off.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from math import gcd

from .formula import molar_mass, parse_formula

__all__ = ["Reaction", "balance", "limiting_reagent"]

_ARROW = re.compile(r"\s*(?:->|=>|→|=)\s*")


def _split_side(side: str) -> list[str]:
    species = [s.strip() for s in re.split(r"\s\+\s|\s*\+\s+|\s+\+\s*", side) if s.strip()]
    out = []
    for s in species:
        m = re.match(r"^(\d+)\s*([A-Z(\[].*)$", s)
        out.append(m.group(2) if m else s)  # ignore any coefficients given by the user
    return out


def _nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    rows, cols = len(matrix), len(matrix[0])
    A = [row[:] for row in matrix]
    pivots = []
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fcol in free:
        v = [Fraction(0)] * cols
        v[fcol] = Fraction(1)
        for i, pc in enumerate(pivots):
            v[pc] = -A[i][fcol]
        basis.append(v)
    return basis


@dataclass(frozen=True)
class Reaction:
    reactants: tuple[str, ...]
    products: tuple[str, ...]
    coefficients: tuple[int, ...]

    @property
    def species(self) -> tuple[str, ...]:
        return self.reactants + self.products

    def coefficient(self, species: str) -> int:
        return self.coefficients[self.species.index(species)]

    def __str__(self) -> str:
        def side(names, coefs):
            return " + ".join(f"{c if c != 1 else ''}{n}" for n, c in zip(names, coefs))
        k = len(self.reactants)
        return f"{side(self.reactants, self.coefficients[:k])} -> {side(self.products, self.coefficients[k:])}"


def balance(equation: str) -> Reaction:
    """Balance an equation such as ``"C3H8 + O2 -> CO2 + H2O"``.

    Raises ``ValueError`` if the equation can't be balanced or has more than one
    independent balanced form.
    """
    sides = _ARROW.split(equation.strip())
    if len(sides) != 2:
        raise ValueError("equation needs exactly one arrow (->, =>, = or →)")
    reactants, products = _split_side(sides[0]), _split_side(sides[1])
    if not reactants or not products:
        raise ValueError("both sides need at least one species")
    species = reactants + products
    comps = [parse_formula(s) for s in species]
    elements = sorted(set().union(*comps))
    sign = [1] * len(reactants) + [-1] * len(products)
    M = [[Fraction(sign[j] * comps[j].get(el, 0)) for j in range(len(species))] for el in elements]

    basis = _nullspace(M)
    if len(basis) == 0:
        raise ValueError("equation cannot be balanced")
    if len(basis) > 1:
        raise ValueError("equation has multiple independent balancings (combine reactions first)")
    v = basis[0]
    lcm = reduce(lambda a, b: a * b // gcd(a, b), (x.denominator for x in v), 1)
    ints = [int(x * lcm) for x in v]
    if ints[0] < 0:
        ints = [-x for x in ints]
    if any(x <= 0 for x in ints):
        raise ValueError("no balancing with all-positive coefficients")
    g = reduce(gcd, ints)
    ints = [x // g for x in ints]
    return Reaction(tuple(reactants), tuple(products), tuple(ints))


def limiting_reagent(reaction: Reaction, masses_g: dict[str, float]) -> dict:
    """Identify the limiting reagent from available reactant masses (grams).

    Returns a dict with the limiting species, the reaction extent (mol), the
    theoretical yield of each product in grams, and the excess remaining of
    each other supplied reactant.
    """
    for s in masses_g:
        if s not in reaction.reactants:
            raise ValueError(f"{s!r} is not a reactant")
    extents = {
        s: (m / molar_mass(s)) / reaction.coefficient(s) for s, m in masses_g.items()
    }
    limiting = min(extents, key=extents.get)
    xi = extents[limiting]
    yields = {p: xi * reaction.coefficient(p) * molar_mass(p) for p in reaction.products}
    excess = {
        s: m - xi * reaction.coefficient(s) * molar_mass(s)
        for s, m in masses_g.items() if s != limiting
    }
    return {"limiting": limiting, "extent_mol": xi, "theoretical_yield_g": yields, "excess_g": excess}
