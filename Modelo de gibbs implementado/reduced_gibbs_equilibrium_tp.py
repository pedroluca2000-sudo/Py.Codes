# coding: utf-8
"""
Reduced-Gibbs Equilibrium (assigned T,P) — π-only solver
---------------------------------------------------------

Implements the reduced Gibbs formulation described in the attached report
(“Implementation of reduced Gibbs iteration equations for chemical equilibrium compositions”,
Invictus II — Propulsion Dept., FEUP, Feb 24, 2026). The core idea is to solve for
Lagrange multipliers π (one per element) and reconstruct the equilibrium composition
via a stabilized softmax. See §7.2–§7.6 of the report for equations and discussion.

Notes
-----
* Standard-state thermodynamics are evaluated with NASA7 polynomials (two ranges),
  using the conventional dimensionless forms h°/(R_u T) and s°/R_u. The standard-state
  reference pressure P° is 1.0e5 Pa by default.
* The solver does not require the total moles n as a Newton unknown; instead, n is
  reconstructed from a reference element r: n = b_r / \\bar{a}_r(π), where \\bar{a} is
  the element-average under the current composition (see §7.4–§7.5 in the report).
* No condensed species: NG = NS (as assumed in the report’s initial implementation).

Usage
-----
1) Prepare a species database (Python dict or YAML) with NASA7 coefficients and
   elemental compositions for the 10 species you want to consider.
2) Build the element totals vector `b` (e.g., from your feed/reactants), covering the
   same element set you pass in `elements` (e.g., ["C","H","O","N"]).
3) Call `solve_equilibrium_tp(T, P, species_db, species_list, elements, b)`.

This module also includes a small CLI demo (`__main__`) that reads a JSON file
with `b` and runs the equilibrium for the default 10-species CHON set.

Author: M365 Copilot for Pedro Laguna (FEUP)
Date: 2026-02-24
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable
import math
import json
import numpy as np

R_u = 8.31446261815324  # J/mol/K (CODATA 2018)

# ---------------------------
# NASA7 thermochemistry
# ---------------------------
@dataclass
class NASA7:
    T_low: float
    T_mid: float
    T_high: float
    coeffs_low: Tuple[float, float, float, float, float, float, float]
    coeffs_high: Tuple[float, float, float, float, float, float, float]

    def region(self, T: float) -> Tuple[float, float, float, float, float, float, float]:
        """Return the 7 coefficients for the appropriate T interval.
        Assumes T_low <= T <= T_high and T_mid splits the low/high ranges.
        """
        if T < self.T_mid:
            return self.coeffs_low
        else:
            return self.coeffs_high

    def cp_R(self, T: float) -> float:
        """Dimensionless cp/R at temperature T (K)."""
        a1, a2, a3, a4, a5, a6, a7 = self.region(T)
        t = T
        return a1 + a2*t + a3*t**2 + a4*t**3 + a5*t**4

    def h_RT(self, T: float) -> float:
        """Dimensionless h°/(R T) at temperature T (K)."""
        a1, a2, a3, a4, a5, a6, a7 = self.region(T)
        t = T
        return a1 + a2*t/2.0 + a3*t**2/3.0 + a4*t**3/4.0 + a5*t**4/5.0 + a6/t

    def s_R(self, T: float) -> float:
        """Dimensionless s°/R at temperature T (K)."""
        a1, a2, a3, a4, a5, a6, a7 = self.region(T)
        t = T
        # Note: a7 is the integration constant for entropy
        return a1*math.log(t) + a2*t + a3*t**2/2.0 + a4*t**3/3.0 + a5*t**4/4.0 + a7

    def g_RT(self, T: float) -> float:
        """Dimensionless g°/(R T) = h°/(R T) - s°/R."""
        return self.h_RT(T) - self.s_R(T)

@dataclass
class Species:
    name: str
    elements: Dict[str, int]  # e.g., {"C":1, "O":2}
    thermo: NASA7

# ---------------------------
# Utilities to build A, g/RT, etc.
# ---------------------------

def build_element_matrix(species_list: List[Species], elements: List[str]) -> np.ndarray:
    """Build A (ℓ x N) where A[i,j] = a_ij (# of atoms of element i in species j)."""
    l = len(elements)
    N = len(species_list)
    A = np.zeros((l, N), dtype=float)
    for j, sp in enumerate(species_list):
        for i, el in enumerate(elements):
            A[i, j] = sp.elements.get(el, 0)
    return A


def g_over_RT_vector(species_list: List[Species], T: float) -> np.ndarray:
    """Vector of g°/(R T) for all species."""
    return np.array([sp.thermo.g_RT(T) for sp in species_list], dtype=float)


def stabilized_softmax(logw: np.ndarray) -> np.ndarray:
    """Compute y = softmax(logw) with overflow protection."""
    m = np.max(logw)
    w = np.exp(logw - m)
    s = np.sum(w)
    return w / s

# ---------------------------
# Reduced-Gibbs π-only solve
# ---------------------------

def solve_equilibrium_tp(
   
    T: float = 2800.0,      # K
    P: float = 3.0e6,       # Pa  (30 bar)
    species_db: Dict[str, Species] = None,
    species_names: List[str] = None,
    elements: List[str] = None,
    b: Dict[str, float] = None,
    P0: float = 1.0e5,      # Pa (pressão padrão)
    max_iter: int = 200,
    tol: float = 1e-10,
    verbose: bool = False,


) -> Dict[str, float]:
    """Solve chemical equilibrium composition at assigned (T,P) for a given element budget.

    Parameters
    ----------
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    species_db : dict[str, Species]
        Database of species objects indexed by name.
    species_names : list[str]
        Species to include in the equilibrium set (e.g., 10 species).
    elements : list[str]
        Element symbols to constrain (e.g., ["C","H","O","N"]).
    b : dict[str, float]
        Element totals (same order as `elements`). Units: kmol-atoms in any consistent
        normalization. Only the ratios matter for composition.
    P0 : float
        Standard-state pressure (Pa) used in the chemical potential. Default 1.0e5 Pa.
    max_iter : int
        Maximum Newton iterations.
    tol : float
        Convergence tolerance on ||f||_2.
    verbose : bool
        If True, prints iteration diagnostics.

    Returns
    -------
    y : dict[name -> mole fraction]
        Equilibrium mole fractions for each species in `species_names`.
    """
    # Assemble working arrays
    sp_list = [species_db[name] for name in species_names]
    A = build_element_matrix(sp_list, elements)  # (ℓ x N)
    g_RT = g_over_RT_vector(sp_list, T)          # (N,)

    ell, N = A.shape
    if ell == 0 or N == 0:
        raise ValueError("Empty element set or species set.")

    b_vec = np.array([b.get(el, 0.0) for el in elements], dtype=float)
    if np.any(b_vec < 0):
        raise ValueError("Element totals 'b' must be nonnegative.")
    if np.allclose(b_vec, 0.0):
        raise ValueError("At least one element total in 'b' must be positive.")

    # Choose reference element r with the largest b_i
    r = int(np.argmax(b_vec))

    # Initialize π with zeros (can be improved with warm starts)
    pi = np.zeros(ell, dtype=float)

    logPterm = -math.log(max(P, 1e-300)/P0)  # appears additively in log-weights

    # Levenberg–Marquardt damping
    lam = 1e-6
    
    def residual_and_jac(pi_vec: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, float]:
        """Compute residual f(π) and analytic Jacobian J(π).
        Returns (f, J, y, n, a_bar, norm_f).
        """
        # log-weights: ℓ_j = -g/RT - ln(P/P0) - π^T a_j
        logw = -g_RT + logPterm - (pi_vec @ A)
        y = stabilized_softmax(logw)  # (N,)
        a_bar = A @ y                 # (ℓ,)
        # total moles from reference element
        a_bar_r = a_bar[r]
        if a_bar_r <= 0:
            # avoid division by zero; regularize slightly
            a_bar_r = 1e-300
        n = b_vec[r] / a_bar_r
        f = n * a_bar - b_vec         # (ℓ,)

        # Build S = A diag(y) A^T efficiently
        # We'll compute A * diag(sqrt(y)) first to improve numerical stability
        sqrt_y = np.sqrt(y)
        A_scaled = A * sqrt_y  # broadcasting column-wise
        S = A_scaled @ A_scaled.T
        # d a_bar / d π_k = -(S_ik - a_bar_i * a_bar_k)
        # Compute Jacobian J_{ik} = dn/dπ_k * a_bar_i + n * d a_bar_i/dπ_k
        # where dn/dπ_k = -(b_r / a_bar_r^2) * d a_bar_r/dπ_k
        d_ab_r_dpi = -(S[:, r] - a_bar * a_bar[r])  # vector over i for element r column
        dn_dpi = -(b_vec[r] / (a_bar_r * a_bar_r)) * d_ab_r_dpi  # (ℓ,)
        # Now full Jacobian
        J = np.empty((ell, ell), dtype=float)
        # Precompute d a_bar / d π as a tensor-free operation using S and outer(a_bar, a_bar)
        # dAbar_dpi[k] = column k of d a_bar / d π = -(S[:,k] - a_bar * a_bar[k])
        for k in range(ell):
            dAbar_dpi_k = -(S[:, k] - a_bar * a_bar[k])  # (ℓ,)
            J[:, k] = dn_dpi[k] * a_bar + n * dAbar_dpi_k
        norm_f = float(np.linalg.norm(f))
        return f, J, y, n, a_bar, norm_f

    f, J, y, n, a_bar, norm_f = residual_and_jac(pi)
    if verbose:
        print(f"iter=0, ||f||={norm_f:.3e}, lam={lam:.1e}")

    for it in range(1, max_iter+1):
        # LM step: (J^T J + lam I) Δπ = -J^T f
        JTJ = J.T @ J
        rhs = -J.T @ f
        # Regularize
        JTJ_reg = JTJ + lam * np.eye(ell)
        try:
            d_pi = np.linalg.solve(JTJ_reg, rhs)
        except np.linalg.LinAlgError:
            # Fallback to least-squares
            d_pi, *_ = np.linalg.lstsq(JTJ_reg, rhs, rcond=None)
        # Backtracking line search
        alpha = 1.0
        success = False
        base_norm = norm_f
        for _bt in range(20):
            pi_trial = pi + alpha * d_pi
            f_trial, J_trial, y_trial, n_trial, a_bar_trial, norm_trial = residual_and_jac(pi_trial)
            if norm_trial < base_norm:
                # accept
                pi = pi_trial
                f, J, y, n, a_bar, norm_f = f_trial, J_trial, y_trial, n_trial, a_bar_trial, norm_trial
                lam = max(lam * 0.3, 1e-12)  # relax damping when successful
                success = True
                break
            else:
                alpha *= 0.5
        if verbose:
            print(f"iter={it}, ||f||={norm_f:.3e}, alpha={alpha:.2e}, lam={lam:.1e}")
        if not success:
            # increase damping and continue
            lam = min(lam * 10.0, 1e12)
        # Convergence check
        if norm_f < tol:
            break
    # Final composition
    y_out = {name: float(val) for name, val in zip(species_names, y)}
    return y_out

# ---------------------------
# Helpers to load a minimal species DB from YAML/JSON or hardcode
# ---------------------------

def make_species_db_from_dict(db: Dict) -> Dict[str, Species]:
    """Build a species_db from a Python dict.

    Expected schema per species:
    {
      "name": "CO2",
      "elements": {"C":1, "O":2},
      "thermo": {
          "T_low": 200.0, "T_mid": 1000.0, "T_high": 6000.0,
          "coeffs_low":  [a1..a7],
          "coeffs_high": [a1..a7]
      }
    }
    """
    species_db: Dict[str, Species] = {}
    for name, item in db.items():
        th = item["thermo"]
        nasa = NASA7(
            T_low=float(th["T_low"]),
            T_mid=float(th["T_mid"]),
            T_high=float(th["T_high"]),
            coeffs_low=tuple(float(x) for x in th["coeffs_low"]),
            coeffs_high=tuple(float(x) for x in th["coeffs_high"]),
        )
        sp = Species(
            name=name,
            elements={k: int(v) for k, v in item["elements"].items()},
            thermo=nasa,
        )
        species_db[name] = sp
    return species_db

# ---------------------------
# __main__ demo (reads JSON files)
# ---------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Reduced-Gibbs (T,P) equilibrium for a fixed element budget (CHON).")
    parser.add_argument("T", type=float, help="Temperature [K]")
    parser.add_argument("P", type=float, help="Pressure [Pa]")
    parser.add_argument("db_json", type=str, help="Path to JSON file with species database (NASA7 + elements)")
    parser.add_argument("b_json", type=str, help="Path to JSON file with element totals, e.g., {\"C\":...,\"H\":...,\"O\":...,\"N\":...}")
    parser.add_argument("--species", nargs="*", default=["CO2","H2O","CO","H2","O2","N2","OH","O","NO","H"],
                        help="Species list (default: 10-species CHON set)")
    parser.add_argument("--elements", nargs="*", default=["C","H","O","N"], help="Element list/order (default: CHON)")
    parser.add_argument("--P0", type=float, default=1.0e5, help="Standard-state pressure [Pa] (default 1e5)")
    parser.add_argument("--tol", type=float, default=1e-10, help="Convergence tolerance on ||f|| (default 1e-10)")
    parser.add_argument("--max_iter", type=int, default=200, help="Max Newton iterations (default 200)")
    parser.add_argument("--verbose", action="store_true", help="Print iteration diagnostics")
    args = parser.parse_args()

    with open(args.db_json, "r", encoding="utf-8") as f:
        db_dict = json.load(f)
    with open(args.b_json, "r", encoding="utf-8") as f:
        b_dict = json.load(f)

    species_db = make_species_db_from_dict(db_dict)
    y = solve_equilibrium_tp(
        T=args.T,
        P=args.P,
        species_db=species_db,
        species_names=args.species,
        elements=args.elements,
        b=b_dict,
        P0=args.P0,
        max_iter=args.max_iter,
        tol=args.tol,
        verbose=args.verbose,
    )
    print(json.dumps(y, indent=2))
