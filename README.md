# stoichiometry-kit

![tests](https://github.com/luisrowork1992-netizen/stoichiometry-kit/actions/workflows/tests.yml/badge.svg)

`stoichkit` is a pure-Python toolkit for everyday **chemical stoichiometry**. It has no dependencies beyond the standard library.

## Features

- **Formula parser** with nested `()`, `[]` and `{}` groups, and hydrates or adducts (`CuSO4·5H2O`, `CaCl2*2H2O`, `K4[Fe(CN)6]`)
- Molar mass and elemental mass fractions, using IUPAC abridged standard atomic weights
- **Exact equation balancing**: rational Gauss–Jordan elimination over the composition matrix, so there is no floating-point guesswork
- Rejects equations that are impossible or that have more than one independent balancing
- **Limiting reagent**, reaction extent, theoretical yields, and leftover excess
- **Empirical formula** from mass percentages, written in Hill order
- Solution-prep helpers: mass for a target molarity, and C₁V₁ = C₂V₂ dilution

## Install

```bash
git clone https://github.com/luisrowork1992-netizen/stoichiometry-kit.git
cd stoichiometry-kit
pip install -e .
```

## Quick start

```python
import stoichkit as sk

sk.molar_mass("CuSO4·5H2O")                       # 249.68 g/mol

rxn = sk.balance("KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2")
print(rxn)   # 2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 8H2O + 5Cl2

res = sk.limiting_reagent(sk.balance("H2 + O2 -> H2O"), {"H2": 4.0, "O2": 32.0})
res["limiting"], res["theoretical_yield_g"]       # ('H2', {'H2O': 35.75...})

sk.empirical_formula({"C": 40.0, "H": 6.71, "O": 53.29})   # 'CH2O'
```

Arrows can be written as `->`, `=>`, `=` or `→`. Any coefficients you type are ignored, and the equation is re-balanced from scratch.

## Notebook

[`notebooks/stoichiometry_walkthrough.ipynb`](notebooks/stoichiometry_walkthrough.ipynb) covers these examples:

- balancing redox and industrial reactions
- thermite as a limiting-reagent problem
- combustion analysis
- CO₂ produced per kilogram for several fuels

## Tests

```bash
pip install pytest
pytest
```

There are 39 tests. They cover:

- parsing, including error handling
- reference molar masses
- balancing seven equations, among them permanganate and phosphorus production
- conservation of atoms
- mass balance in limiting-reagent results
- empirical formulas

## License

MIT
