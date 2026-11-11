import pytest

import stoichkit as sk


@pytest.mark.parametrize(
    "formula, expected",
    [
        ("H2O", {"H": 2, "O": 1}),
        ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}),
        ("Al2(SO4)3", {"Al": 2, "S": 3, "O": 12}),
        ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}),
        ("CuSO4·5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}),
        ("CaCl2*2H2O", {"Ca": 1, "Cl": 2, "H": 4, "O": 2}),
        ("(CH3)3COH", {"C": 4, "H": 10, "O": 1}),
    ],
)
def test_parse_formula(formula, expected):
    assert dict(sk.parse_formula(formula)) == expected


@pytest.mark.parametrize("bad", ["", "Xx2", "Ca(OH", "H2O)", "2", "(]"])
def test_parse_errors(bad):
    with pytest.raises(ValueError):
        sk.parse_formula(bad)


@pytest.mark.parametrize(
    "formula, mass",
    [
        ("H2O", 18.015),
        ("CO2", 44.009),
        ("NaCl", 58.44),
        ("C6H12O6", 180.156),
        ("CuSO4·5H2O", 249.677),
        ("K4[Fe(CN)6]", 368.345),
    ],
)
def test_molar_mass(formula, mass):
    assert sk.molar_mass(formula) == pytest.approx(mass, abs=0.01)


def test_mass_fractions_sum_to_one():
    fr = sk.mass_fractions("C6H12O6")
    assert sum(fr.values()) == pytest.approx(1.0)
    assert fr["C"] == pytest.approx(0.4000, abs=1e-3)


@pytest.mark.parametrize(
    "eq, coeffs",
    [
        ("H2 + O2 -> H2O", (2, 1, 2)),
        ("C3H8 + O2 -> CO2 + H2O", (1, 5, 3, 4)),
        ("Fe2O3 + CO -> Fe + CO2", (1, 3, 2, 3)),
        ("KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2", (2, 16, 2, 2, 8, 5)),
        ("Al + O2 = Al2O3", (4, 3, 2)),
        ("C8H18 + O2 → CO2 + H2O", (2, 25, 16, 18)),
        ("Ca3(PO4)2 + SiO2 + C -> CaSiO3 + P4 + CO", (2, 6, 10, 6, 1, 10)),
    ],
)
def test_balance(eq, coeffs):
    assert sk.balance(eq).coefficients == coeffs


def test_balanced_reaction_conserves_atoms():
    rxn = sk.balance("KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2")
    k = len(rxn.reactants)
    left, right = {}, {}
    for sp, c in zip(rxn.species, rxn.coefficients):
        target = left if rxn.species.index(sp) < k else right
        for el, n in sk.parse_formula(sp).items():
            target[el] = target.get(el, 0) + c * n
    assert left == right


def test_balance_ignores_user_coefficients():
    assert sk.balance("2H2 + 1O2 -> 2H2O").coefficients == (2, 1, 2)


def test_balance_errors():
    with pytest.raises(ValueError):
        sk.balance("H2 + O2")                      # no arrow
    with pytest.raises(ValueError):
        sk.balance("H2 -> O2")                     # impossible
    with pytest.raises(ValueError):
        sk.balance("H2 + O2 -> H2O + H2O2")         # two independent reactions


def test_reaction_str():
    assert str(sk.balance("H2 + O2 -> H2O")) == "2H2 + O2 -> 2H2O"


def test_limiting_reagent_water():
    rxn = sk.balance("H2 + O2 -> H2O")
    res = sk.limiting_reagent(rxn, {"H2": 4.0, "O2": 32.0})
    n_h2 = 4.0 / sk.molar_mass("H2")
    assert res["limiting"] == "H2"
    assert res["theoretical_yield_g"]["H2O"] == pytest.approx(n_h2 * sk.molar_mass("H2O"))
    assert res["excess_g"]["O2"] == pytest.approx(32.0 - n_h2 / 2 * sk.molar_mass("O2"))


def test_limiting_reagent_mass_balance():
    rxn = sk.balance("Fe2O3 + CO -> Fe + CO2")
    supplied = {"Fe2O3": 100.0, "CO": 30.0}
    res = sk.limiting_reagent(rxn, supplied)
    consumed = sum(supplied.values()) - sum(res["excess_g"].values())
    produced = sum(res["theoretical_yield_g"].values())
    assert consumed == pytest.approx(produced, rel=1e-9)


@pytest.mark.parametrize(
    "pct, formula",
    [
        ({"C": 40.0, "H": 6.71, "O": 53.29}, "CH2O"),
        ({"C": 85.63, "H": 14.37}, "CH2"),
        ({"Fe": 69.94, "O": 30.06}, "Fe2O3"),
        ({"C": 92.26, "H": 7.74}, "CH"),
        ({"K": 26.58, "Cr": 35.35, "O": 38.07}, "Cr2K2O7"),
    ],
)
def test_empirical_formula(pct, formula):
    assert sk.empirical_formula(pct) == formula


def test_solution_helpers():
    # 0.5 L of 0.1 M NaCl needs 2.922 g
    assert sk.mass_for_molarity("NaCl", 0.1, 0.5) == pytest.approx(2.922, abs=1e-3)
    assert sk.molarity_from_mass("NaCl", 2.922, 0.5) == pytest.approx(0.1, rel=1e-3)
    assert sk.dilution_volume(12.0, 1.0, 250.0) == pytest.approx(20.8333, rel=1e-4)
    with pytest.raises(ValueError):
        sk.dilution_volume(1.0, 2.0, 100.0)
