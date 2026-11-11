"""Solution-preparation helpers."""

from __future__ import annotations

from .formula import molar_mass

__all__ = ["mass_for_molarity", "dilution_volume", "molarity_from_mass"]


def mass_for_molarity(formula: str, molarity: float, volume_L: float) -> float:
    """Grams of solute needed for ``volume_L`` litres of ``molarity`` mol/L solution."""
    return molarity * volume_L * molar_mass(formula)


def molarity_from_mass(formula: str, mass_g: float, volume_L: float) -> float:
    """Molarity (mol/L) of ``mass_g`` grams dissolved to ``volume_L`` litres."""
    return mass_g / molar_mass(formula) / volume_L


def dilution_volume(c_stock: float, c_target: float, v_target: float) -> float:
    """Stock volume needed for a dilution, from ``C1 V1 = C2 V2`` (same volume unit as ``v_target``)."""
    if c_target > c_stock:
        raise ValueError("target concentration exceeds the stock concentration")
    return c_target * v_target / c_stock
