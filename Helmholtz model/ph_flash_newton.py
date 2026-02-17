# ph_flash_newton.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Dict, Tuple
import math

class EOSBackendProto:
    def pressure(self, rho_mass: float, T: float) -> float: ...
    def enthalpy_mass(self, rho_mass: float, T: float) -> float: ...

@dataclass
class PHFlashOptions:
    tol_f: float = 1e-10
    tol_x: float = 1e-12
    max_iter: int = 50
    line_search: bool = True
    ls_beta: float = 0.5
    ls_max: int = 12
    enforce_physical: bool = True
    fd_rel_step_rho: float = 1e-4
    fd_rel_step_T: float = 1e-4

def _num_jac(backend: EOSBackendProto, rho: float, T: float, Pstar: float, Hstar: float, opts: PHFlashOptions):
    def f_eval(r,t):
        return (backend.pressure(r,t) - Pstar, backend.enthalpy_mass(r,t) - Hstar)
    fr, hr = f_eval(rho, T)

    dr = max(abs(rho),1.0)*opts.fd_rel_step_rho
    dT = max(abs(T),  1.0)*opts.fd_rel_step_T
    if opts.enforce_physical:
        if rho - dr <= 0: dr = 0.5*rho
        if T   - dT <= 0: dT = 0.5*T

    f_p = f_eval(rho+dr, T);      f_m = f_eval(rho-dr, T)
    dp_dr = (f_p[0]-f_m[0])/(2*dr);  dh_dr = (f_p[1]-f_m[1])/(2*dr)

    f_p = f_eval(rho, T+dT);      f_m = f_eval(rho, T-dT)
    dp_dT = (f_p[0]-f_m[0])/(2*dT);  dh_dT = (f_p[1]-f_m[1])/(2*dT)

    J = ((dp_dr, dp_dT), (dh_dr, dh_dT))
    f = (fr, hr)
    return J, f

def _solve_2x2(J, f):
    (a,b), (c,d) = J
    det = a*d - b*c
    if not (det != 0.0):  # cobre 0 e NaN
        raise RuntimeError("Jacobiano quase singular.")
    inv = ((d/det, -b/det), (-c/det, a/det))
    dx0 = -(inv[0][0]*f[0] + inv[0][1]*f[1])
    dx1 = -(inv[1][0]*f[0] + inv[1][1]*f[1])
    return dx0, dx1

def ph_flash(backend: EOSBackendProto, P_star: float, H_star: float,
             T0: float, rho0: float, options: Optional[PHFlashOptions]=None) -> Tuple[float,float,Dict]:
    opts = options or PHFlashOptions()
    rho = max(rho0, 1e-10); T = max(T0, 1e-10)
    hist = []

    for it in range(opts.max_iter):
        J, f = _num_jac(backend, rho, T, P_star, H_star, opts)
        fp_rel = abs(f[0])/max(abs(P_star),1.0)
        fh_rel = abs(f[1])/max(abs(H_star),1.0)
        hist.append({"iter":it,"rho":rho,"T":T,"fp":f[0],"fh":f[1],"frel":max(fp_rel, fh_rel)})

        if max(fp_rel, fh_rel) < opts.tol_f:
            return rho, T, {"converged":True, "it":it, "history":hist}

        drho, dT = _solve_2x2(J, f)

        small = max(abs(drho)/max(abs(rho),1.0), abs(dT)/max(abs(T),1.0)) < opts.tol_x
        if small:
            return rho, T, {"converged":True, "it":it, "history":hist}

        step = 1.0; accepted = False
        for _ in range(opts.ls_max if opts.line_search else 1):
            rr = rho + step*drho; TT = T + step*dT
            if opts.enforce_physical and (rr<=0.0 or TT<=0.0):
                step *= opts.ls_beta; continue
            fp = backend.pressure(rr, TT) - P_star
            fh = backend.enthalpy_mass(rr, TT) - H_star
            frel_try = max(abs(fp)/max(abs(P_star),1.0), abs(fh)/max(abs(H_star),1.0))
            if frel_try < max(fp_rel, fh_rel):
                rho, T = rr, TT; accepted = True; break
            step *= opts.ls_beta

        if not accepted:
            rho = max(rho + step*drho, 1e-20) if opts.enforce_physical else rho + step*drho
            T   = max(T   + step*dT,   1e-20) if opts.enforce_physical else T   + step*dT

    return rho, T, {"converged":False, "it":opts.max_iter, "history":hist}