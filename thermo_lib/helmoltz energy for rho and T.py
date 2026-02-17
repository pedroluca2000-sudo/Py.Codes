"""
Flash P-H (encontrar T e rho dado P*, H*) com Newton 2x2.
Base mássica (SI): P [Pa], T [K], rho [kg/m3], h [J/kg].

Arquitetura:
- EOSBackend: interface mínima (p(rho,T), h(rho,T)).
- IdealGasBackend: backend de demonstração (gás ideal com c_p constante).
- ph_flash(): resolvedor genérico com Jacobiano numérico (dif. centrais).

Para plugar uma EoS de Helmholtz (IAPWS-95, Span–Wagner), implemente um
backend que calcule p(rho,T) e h(rho,T) a partir de α(δ,τ) e suas derivadas:
    p = rho*R*T*(1 + δ*α_δ^r)
    h = R*T*( 1 + δ*α_δ^r + τ*(α_τ^0 + α_τ^r) )   [em base molar]
e converta para base mássica dividindo por M. Ver:
- IAPWS-95 (água) e suas fórmulas de derivadas. [Wagner & Pruß, 2002]
- Span–Wagner (CO2). [Span & Wagner, 1996]
- Resumo de derivadas (teqp, NIST). [teqp docs]
"""

from dataclasses import dataclass
from typing import Tuple, Dict, Optional, Callable
import math

# ============================================================
# 1) Backend termodinâmico: interface mínima
# ============================================================

class EOSBackend:
    """Interface mínima: fornece p(rho,T) e h(rho,T) em base mássica (SI)."""
    def p(self, rho: float, T: float) -> float:
        """Pressão [Pa]"""
        raise NotImplementedError

    def h(self, rho: float, T: float) -> float:
        """Entalpia específica [J/kg]"""
        raise NotImplementedError


# ------------------------------------------------------------
# 1a) Backend de demonstração: GÁS IDEAL (c_p constante)
#     p = rho * R * T
#     h(T) = h_ref + c_p * (T - T_ref)
# ------------------------------------------------------------

@dataclass
class IdealGasBackend(EOSBackend):
    R: float = 287.0        # [J/(kg.K)] -> ar seco ~ 287
    cp: float = 1006.0      # [J/(kg.K)] -> ar seco ~ 1006
    T_ref: float = 298.15   # [K]
    h_ref: float = 0.0      # [J/kg] entalpia na referência (ajustável)

    def p(self, rho: float, T: float) -> float:
        return rho * self.R * T

    def h(self, rho: float, T: float) -> float:
        return self.h_ref + self.cp * (T - self.T_ref)


# ============================================================
# 2) Resolvedor P-H (Newton 2x2 com Jacobiano numérico)
# ============================================================

tole = 

@dataclass
class PHFlashOptions:
    tol_f: float = 1e-8        # tolerância relativa nos resíduos
    tol_x: float = 1e-10       # tolerância relativa nos incrementos
    max_iter: int = 50
    fd_rel_step_rho: float = 1e-6
    fd_rel_step_T: float = 1e-6
    line_search: bool = True
    ls_beta: float = 0.5       # fator de recuo da busca linear
    ls_max: int = 8            # máx. passos na busca linear
    enforce_physical: bool = True  # manter T>0, rho>0


def _numerical_jacobian(
    backend: EOSBackend,
    rho: float, T: float,
    P_star: float, H_star: float,
    opts: PHFlashOptions
) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """
    Jacobiano numérico J = [[dp/drho, dp/dT], [dh/drho, dh/dT]]
    e vetores f = [p-P*, h-H*] e (p,h) atuais.
    Diferenças centrais com passos relativos em rho e T.
    """
    # Funções
    def f_eval(r, t):
        p = backend.p(r, t)
        h = backend.h(r, t)
        return (p - P_star, h - H_star), p, h

    # Passos de FD
    drho = max(abs(rho), 1.0) * opts.fd_rel_step_rho
    dT   = max(abs(T),   1.0) * opts.fd_rel_step_T

    # Ajustes para manter positividade se pedido
    if opts.enforce_physical:
        if rho - drho <= 0.0: drho = 0.5 * rho
        if T   - dT   <= 0.0: dT   = 0.5 * T

    # Estado central
    f0, p0, h0 = f_eval(rho, T)

    # Derivadas em rho (T constante)
    f_plus, _, _  = f_eval(rho + drho, T)
    f_minus, _, _ = f_eval(rho - drho, T)
    dp_drho = (f_plus[0] - f_minus[0]) / (2.0 * drho)
    dh_drho = (f_plus[1] - f_minus[1]) / (2.0 * drho)

    # Derivadas em T (rho constante)
    f_plus, _, _  = f_eval(rho, T + dT)
    f_minus, _, _ = f_eval(rho, T - dT)
    dp_dT = (f_plus[0] - f_minus[0]) / (2.0 * dT)
    dh_dT = (f_plus[1] - f_minus[1]) / (2.0 * dT)

    J = ((dp_drho, dp_dT),
         (dh_drho, dh_dT))

    return J, f0, (p0, h0)


def _solve_2x2(J, f):
    """Resolve J * dx = -f para dx = [drho, dT]."""
    (a, b), (c, d) = J
    det = a * d - b * c
    if abs(det) < 1e-30:
        raise RuntimeError("Jacobiano quase singular.")
    invJ = (( d / det, -b / det),
            (-c / det,  a / det))
    # dx = -invJ @ f
    drho = -(invJ[0][0] * f[0] + invJ[0][1] * f[1])
    dT   = -(invJ[1][0] * f[0] + invJ[1][1] * f[1])
    return drho, dT


def ph_flash(
    backend: EOSBackend,
    P_star: float, H_star: float,
    T0: float, rho0: float,
    options: Optional[PHFlashOptions] = None
) -> Tuple[float, float, Dict]:
    """
    Resolve p(rho,T)=P*, h(rho,T)=H* para (T, rho).
    Retorna (rho, T, diag).
    """
    opts = options or PHFlashOptions()

    # Inicialização
    rho = max(rho0, 1e-12)
    T   = max(T0,   1e-12)

    history = []
    for it in range(opts.max_iter):
        J, f, (p_curr, h_curr) = _numerical_jacobian(backend, rho, T, P_star, H_star, opts)
        # Critérios de convergência (resíduos relativos)
        f1_rel = abs(f[0]) / max(abs(P_star), 1.0)
        f2_rel = abs(f[1]) / max(abs(H_star), 1.0)
        history.append({
            "iter": it, "rho": rho, "T": T,
            "p": p_curr, "h": h_curr,
            "f_p": f[0], "f_h": f[1],
            "f_rel": max(f1_rel, f2_rel)
        })
        if max(f1_rel, f2_rel) < opts.tol_f:
            return rho, T, {"converged": True, "it": it, "history": history}

        # Passo de Newton
        drho, dT = _solve_2x2(J, f)

        # Critério de passo pequeno
        if max(abs(drho) / max(abs(rho), 1.0),
               abs(dT)   / max(abs(T),   1.0)) < opts.tol_x:
            return rho, T, {"converged": True, "it": it, "history": history}

        # Busca linear (opcional) para robustez
        step = 1.0
        accepted = False
        for _ in range(opts.ls_max if opts.line_search else 1):
            rho_try = rho + step * drho
            T_try   = T   + step * dT
            if opts.enforce_physical:
                if rho_try <= 0.0 or T_try <= 0.0:
                    step *= opts.ls_beta
                    continue
            # Avaliar novo resíduo
            f_try_p = backend.p(rho_try, T_try) - P_star
            f_try_h = backend.h(rho_try, T_try) - H_star
            f_try_rel = max(abs(f_try_p)/max(abs(P_star),1.0),
                            abs(f_try_h)/max(abs(H_star),1.0))
            if f_try_rel < max(f1_rel, f2_rel):
                rho, T = rho_try, T_try
                accepted = True
                break
            step *= opts.ls_beta

        if not accepted:
            # Se nenhuma melhora, ainda assim dá o passo (pequeno)
            rho = max(rho + step * drho, 1e-20) if opts.enforce_physical else rho + step * drho
            T   = max(T   + step * dT,   1e-20) if opts.enforce_physical else T   + step * dT

    return rho, T, {"converged": False, "it": opts.max_iter, "history": history}


# ============================================================
# 3) Exemplo de uso com GÁS IDEAL (teste de sanidade)
# ============================================================

if __name__ == "__main__":
    # "Ar seco" aproximado:
    gas = IdealGasBackend(R=287.0, cp=1006.0, T_ref=298.15, h_ref=0.0)

    # Estado verdadeiro:
    P_true = 1.0e5             # Pa
    T_true = 300.0             # K
    rho_true = P_true / (gas.R * T_true)
    H_true = gas.h(rho_true, T_true)

    # Chutes:
    T0 = 320.0
    rho0 = 1.2  # kg/m3 (ordem de grandeza do ar)

    rho_sol, T_sol, diag = ph_flash(gas, P_true, H_true, T0, rho0)

    print("Convergência:", diag["converged"], "iterações:", diag["it"])
    print(f"Resultado: T = {T_sol:.6f} K, rho = {rho_sol:.8f} kg/m3")
    print(f"Verdadeiro: T = {T_true:.6f} K, rho = {rho_true:.8f} kg/m3")