import json
import os
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = Path(os.environ.get(
    "OGREPY_NOTEBOOK", ROOT / "Einstein_Euler_BDNK_Calculator.ipynb"
))

os.environ["OGREPY_DISABLE_WELCOME"] = "True"


def run_selected_cells():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    namespace = {"__name__": "__main__"}
    # The full conservation-law and Einstein-BDNK expansions in cell 21 are
    # deliberately excluded from the fast algebraic verification pass. They
    # are very large, while the constitutive identities below do not depend on
    # expanding those PDEs component-by-component.
    selected = [1, 3, 4, 6, 7, 9, 13, 15, 18, 23]
    for index in selected:
        source = "".join(notebook["cells"][index].get("source", []))
        start = time.perf_counter()
        print(f"RUN cell {index}", flush=True)
        exec(compile(source, f"cell_{index}", "exec"), namespace)
        elapsed = time.perf_counter() - start
        print(f"DONE cell {index}: {elapsed:.2f}s", flush=True)
    return namespace


def tensor_all_zero(ns, tensor, indices):
    ok, values = ns["all_zero_components"](tensor, indices=indices, simplify=True)
    return bool(ok), values


def sympy_array_all_zero(ns, array):
    og = ns["og"]
    values = ns["flatten_array"](array)
    simplified = [og.s.simplify(value) for value in values]
    return all(value == 0 for value in simplified), simplified


def run_independent_checks(ns):
    og = ns["og"]
    mu = ns["mu"]
    nu = ns["nu"]
    alpha = og.sym("alpha")
    beta = og.sym("beta")
    EF = ns["EF"]

    checks = {}

    velocity_norm = og.calc(
        formula=ns["Velocity"](mu) @ ns["Velocity"](mu),
        symbol=r"u\cdot u",
    )
    norm_value = og.s.simplify(
        velocity_norm.components(coords=EF, indices=(), warn=False)[0] + 1
    )
    checks["u_mu u^mu = -1"] = norm_value == 0

    q_orth = og.calc(
        formula=ns["Velocity"](mu) @ ns["Q"](mu),
        symbol=r"u\cdot Q",
    )
    q_value = og.s.simplify(q_orth.components(coords=EF, indices=(), warn=False)[0])
    checks["u_mu Q^mu = 0"] = q_value == 0

    j_orth = og.calc(
        formula=ns["Velocity"](mu) @ ns["Jdiss"](mu),
        symbol=r"u\cdot J_{diss}",
    )
    j_value = og.s.simplify(j_orth.components(coords=EF, indices=(), warn=False)[0])
    checks["u_mu Jdiss^mu = 0"] = j_value == 0

    number_projection = og.calc(
        formula=ns["Velocity"](mu) @ ns["BDNKCurrent"](mu),
        symbol=r"u\cdot J_{BDNK}",
    )
    number_value = og.s.simplify(
        number_projection.components(coords=EF, indices=(), warn=False)[0]
        + ns["Ncal_expr"]
    )
    checks["-u_mu J_BDNK^mu = Ncal"] = number_value == 0

    current_spatial_projection = og.calc(
        formula=(ns["Delta"](mu, alpha) @ ns["BDNKCurrent"](alpha)) - ns["Jdiss"](mu),
        symbol=r"Delta J-J_{diss}",
    )
    checks["Delta^mu_nu J_BDNK^nu = Jdiss^mu"] = tensor_all_zero(
        ns, current_spatial_projection, indices=(1,)
    )[0]

    stress_sym = og.calc(
        formula=ns["BDNKStress"](mu, nu) - ns["BDNKStress"](nu, mu),
        symbol=r"T-T^T",
    )
    checks["T_BDNK symmetry"] = tensor_all_zero(ns, stress_sym, indices=(1, 1))[0]

    energy_projection = og.calc(
        formula=ns["Velocity"](mu) @ ns["BDNKStress"](mu, nu) @ ns["Velocity"](nu),
        symbol=r"uTu",
    )
    energy_value = og.s.simplify(
        energy_projection.components(coords=EF, indices=(), warn=False)[0]
        - ns["Ecal_expr"]
    )
    checks["u_mu u_nu T_BDNK^{mu nu} = Ecal"] = energy_value == 0

    pressure_projection = og.calc(
        formula=ns["Delta"](mu, nu) @ ns["BDNKStress"](mu, nu),
        symbol=r"Delta T",
    )
    pressure_value = og.s.simplify(
        pressure_projection.components(coords=EF, indices=(), warn=False)[0]
        - 3 * ns["Pcal_expr"]
    )
    checks["(1/3) Delta_mn T_BDNK^{mn} = Pcal"] = pressure_value == 0

    heat_flux_projection = og.calc(
        formula=(
            ns["Delta"](mu, alpha)
            @ ns["Velocity"](beta)
            @ ns["BDNKStress"](alpha, beta)
        ) + ns["Q"](mu),
        symbol=r"Delta uT+Q",
    )
    checks["-Delta^mu_a u_b T_BDNK^{ab} = Q^mu"] = tensor_all_zero(
        ns, heat_flux_projection, indices=(1,)
    )[0]

    # In the zero-gradient/zero-viscosity limit, the BDNK constitutive tensors
    # must reduce algebraically to the perfect-fluid tensors built from the
    # same equation of state.
    coeffs = [
        ns[name]
        for name in (
            "eps1", "eps2", "eps3",
            "pi1", "pi2", "pi3",
            "vartheta1", "vartheta2", "vartheta3",
            "nu1", "nu2", "nu3",
            "gamma1", "gamma2", "gamma3",
            "eta",
        )
    ]
    substitutions = {coefficient: 0 for coefficient in coeffs}

    bdnk_stress = ns["BDNKStress"].components(coords=EF, indices=(1, 1), warn=False)
    ideal_stress = og.calc(
        formula=(
            ns["eps0"] * (ns["Velocity"](mu) @ ns["Velocity"](nu))
            + ns["P0"] * ns["Delta"](mu, nu)
        ),
        symbol=r"T_{ideal}",
    ).components(coords=EF, indices=(1, 1), warn=False)
    stress_difference = bdnk_stress.applyfunc(lambda value: value.subs(substitutions)) - ideal_stress
    checks["perfect-fluid stress limit"] = sympy_array_all_zero(ns, stress_difference)[0]

    bdnk_current = ns["BDNKCurrent"].components(coords=EF, indices=(1,), warn=False)
    ideal_current = og.calc(
        formula=ns["n0"] * ns["Velocity"](mu),
        symbol=r"J_{ideal}",
    ).components(coords=EF, indices=(1,), warn=False)
    current_difference = bdnk_current.applyfunc(lambda value: value.subs(substitutions)) - ideal_current
    checks["perfect-fluid current limit"] = sympy_array_all_zero(ns, current_difference)[0]

    stress_components = bdnk_stress
    angular_relation = og.s.simplify(
        stress_components[3, 3] * og.s.sin(ns["theta"]) ** 2 - stress_components[2, 2]
    )
    checks["spherical angular relation T^33 sin^2(theta)=T^22"] = angular_relation == 0

    print("\nINDEPENDENT CHECKS")
    for name, passed in checks.items():
        print(f"{name}: {passed}")

    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise AssertionError(f"Failed checks: {failed}")


if __name__ == "__main__":
    state = run_selected_cells()
    run_independent_checks(state)
