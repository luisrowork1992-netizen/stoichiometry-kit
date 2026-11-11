"""Chemical formula parsing, molar masses and composition."""

from __future__ import annotations

import re
from collections import Counter
from functools import reduce
from math import gcd

from .elements import ATOMIC_WEIGHTS

__all__ = [
    "parse_formula",
    "molar_mass",
    "mass_fractions",
    "empirical_formula",
    "format_formula",
]

_TOKEN = re.compile(r"([A-Z][a-z]?)|(\d+)|([(\[{])|([)\]}])|([·*.])")
_CLOSE = {")": "(", "]": "[", "}": "{"}


def parse_formula(formula: str) -> Counter:
    """Parse a formula into element counts.

    Supports nested groups ``()[]{}``, and hydrates or adducts joined by
    ``·``, ``*`` or ``.`` with an optional leading multiplier, e.g.
    ``CuSO4·5H2O`` or ``K4[Fe(CN)6]``.
    """
    formula = formula.replace(" ", "")
    if not formula:
        raise ValueError("empty formula")
    parts = re.split(r"[·*.]", formula)
    total: Counter = Counter()
    for part in parts:
        m = re.match(r"^(\d+)(.*)$", part)
        mult, body = (int(m.group(1)), m.group(2)) if m else (1, part)
        if not body:
            raise ValueError(f"empty component in {formula!r}")
        for el, n in _parse_group(body).items():
            total[el] += n * mult
    return total


def _parse_group(s: str) -> Counter:
    stack: list[tuple[str, Counter]] = [("", Counter())]
    pos = 0
    last: Counter | None = None  # what the next number multiplies
    while pos < len(s):
        m = _TOKEN.match(s, pos)
        if not m:
            raise ValueError(f"unexpected character {s[pos]!r} in {s!r}")
        el, num, opn, cls, _ = m.groups()
        if el:
            if el not in ATOMIC_WEIGHTS:
                raise ValueError(f"unknown element {el!r}")
            last = Counter({el: 1})
            stack[-1][1].update(last)
        elif num:
            if last is None:
                raise ValueError(f"number without element in {s!r}")
            for k, v in last.items():
                stack[-1][1][k] += v * (int(num) - 1)
            last = None
        elif opn:
            stack.append((opn, Counter()))
            last = None
        elif cls:
            if len(stack) == 1 or stack[-1][0] != _CLOSE[cls]:
                raise ValueError(f"unbalanced brackets in {s!r}")
            _, group = stack.pop()
            stack[-1][1].update(group)
            last = group
        pos = m.end()
    if len(stack) != 1:
        raise ValueError(f"unbalanced brackets in {s!r}")
    return +stack[0][1]


def molar_mass(formula: str) -> float:
    """Molar mass in g/mol."""
    return sum(ATOMIC_WEIGHTS[el] * n for el, n in parse_formula(formula).items())


def mass_fractions(formula: str) -> dict[str, float]:
    """Mass fraction of each element (sums to 1)."""
    counts = parse_formula(formula)
    M = molar_mass(formula)
    return {el: ATOMIC_WEIGHTS[el] * n / M for el, n in counts.items()}


def format_formula(counts: dict[str, int]) -> str:
    """Hill-order formula string (C first, then H, then alphabetical)."""
    keys = sorted(counts)
    if "C" in counts:
        keys = ["C"] + (["H"] if "H" in counts else []) + [k for k in keys if k not in ("C", "H")]
    return "".join(f"{k}{counts[k] if counts[k] != 1 else ''}" for k in keys)


def empirical_formula(mass_percent: dict[str, float], max_multiplier: int = 12, tol: float = 0.05) -> str:
    """Empirical formula from elemental mass percentages.

    Mole ratios are scaled by the smallest integer (up to ``max_multiplier``)
    that brings every ratio within ``tol`` of a whole number.
    """
    moles = {el: pct / ATOMIC_WEIGHTS[el] for el, pct in mass_percent.items()}
    smallest = min(moles.values())
    ratios = {el: n / smallest for el, n in moles.items()}
    for k in range(1, max_multiplier + 1):
        scaled = {el: r * k for el, r in ratios.items()}
        if all(abs(v - round(v)) <= tol * k for v in scaled.values()):
            ints = {el: int(round(v)) for el, v in scaled.items()}
            g = reduce(gcd, ints.values())
            return format_formula({el: v // g for el, v in ints.items()})
    raise ValueError("could not find a small whole-number ratio")
