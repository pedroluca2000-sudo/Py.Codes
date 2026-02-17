# run_ph_flash_N2O.py
from __future__ import annotations
from ph_flash_newton import ph_flash, PHFlashOptions
from params_N2O_lemmon_span_2006 import make_eos

if __name__ == "__main__":
    eos = make_eos()

    # Dados que você passou (P em Pa, h em kJ/kg -> converter):
    P = 4_995_518.404026436
    H = 191_643.86209589397 * 1e3  # kJ/kg -> J/kg

    # Chutes (próximos do estado líquido comprimido ~285 K, ρ ~ 851 kg/m3)
    T0   = 300.0
    rho0 = 800.0

    opts = PHFlashOptions(tol_f=1e-10, tol_x=1e-12, line_search=True)

    rho_sol, T_sol, diag = ph_flash(eos, P, H, T0, rho0, options=opts)

    print("Convergência:", diag["converged"], "iterações:", diag["it"])
    print(f"Solução: T = {T_sol:.9f} K, ρ = {rho_sol:.12f} kg/m3")
    print(f"Check: p = {eos.pressure(rho_sol, T_sol):.6f} Pa, h = {eos.enthalpy_mass(rho_sol, T_sol):.6f} J/kg")