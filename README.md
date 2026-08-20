# Symbolic Einstein–Euler / BDNK Calculator

This project uses [OGRePy](https://github.com/bshoshany/OGRePy) to construct
the Einstein–Euler equations and the general first-order BDNK constitutive
tensors in a spherically symmetric ingoing-null coordinate system.

The metric and four-velocity ansatz are

$$
ds^2=-a(t,r)b(t,r)^2dt^2+2b(t,r)\,dt\,dr+r^2d\Omega^2,
$$

$$
u^\mu=\left(U,\frac{a b^2U^2-1}{2bU},0,0\right),
\qquad g_{\mu\nu}u^\mu u^\nu=-1.
$$

The notebook constructs the ideal-fluid current and stress tensor, the BDNK
geometric building blocks, the general first-order current and stress tensor,
their conservation equations, and the corresponding Einstein equations.

## What has been checked

The constitutive definitions were compared with Eqs. (4)–(6) of
Bemfica–Disconzi–Noronha. Independent symbolic checks verify:

- four-velocity normalization;
- projector and acceleration orthogonality;
- symmetry, tracelessness, and transversality of the shear tensor;
- transversality of $Q^\mu$ and $\mathcal J^\mu$;
- recovery of $\mathcal E$, $\mathcal P$, $\mathcal N$, $Q^\mu$, and
  $\mathcal J^\mu$ by the standard projections of $T^{\mu\nu}$ and $J^\mu$;
- reduction to the perfect-fluid tensors when first-order coefficients vanish;
- the expected angular-component relation from spherical symmetry.

See [`VERIFICATION.md`](VERIFICATION.md) for the exact boundary between what
was verified and what remains unchecked.

## Run

Tested with Python 3.12, OGRePy 1.3.1, and SymPy 1.14.0.

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
jupyter lab Einstein_Euler_BDNK_Calculator.ipynb
```

The full BDNK conservation-law and Einstein-equation expansions are expensive.
They are disabled by default. Set `BUILD_FULL_BDNK_EQUATIONS = True` near the
top of the notebook to construct and export them. Pre-generated LaTeX
components from the checked run are included in [`outputs/`](outputs/).

## Important limitation

This repository implements the **general first-order constitutive ansatz**.
Arbitrary symbolic transport coefficients are not automatically causal,
stable, thermodynamically consistent, or strongly hyperbolic. In particular,
thermodynamic consistency imposes $\vartheta_1=\vartheta_2$ and
$\gamma_1=\gamma_2$ in the notation of the notebook, while causal and stable
frames require additional equation-of-state and coefficient inequalities.

This is a symbolic research calculator, not a numerical BDNK evolution code.

## References

- F. S. Bemfica, M. M. Disconzi, and J. Noronha, *First-Order
  General-Relativistic Viscous Fluid Dynamics*, arXiv:2009.11388.
- B. Shoshany, *OGRePy: An Object-Oriented General Relativity Package for
  Python*, Journal of Open Research Software 13, 9 (2025).
- L. S. Keeble and F. Pretorius, *First-Order Viscous Relativistic
  Hydrodynamics on the Two-Sphere*, Physical Review D 112, 124034 (2025).

Machine-readable entries are in [`references.bib`](references.bib).
