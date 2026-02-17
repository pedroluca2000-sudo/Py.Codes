# eos_helmholtz_core.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, Tuple

# ==========================
# Definições e utilidades
# ==========================

R_u = 8.31446261815324  # J/(mol.K) - universal

class ResidualTerm(Protocol):
    """
    Protocolo de um termo residual da forma short (polinomial/exponencial).
    Cada termo sabe avaliar: φr, φr_δ, φr_τ, φr_δδ, φr_δτ, φr_ττ.
    """
    def val(self, delta: float, tau: float) -> float: ...
    def d_delta(self, delta: float, tau: float) -> float: ...
    def d_tau(self,   delta: float, tau: float) -> float: ...
    def d_delta2(self,delta: float, tau: float) -> float: ...
    def d_deltatau(self,delta: float, tau: float) -> float: ...
    def d_tau2(self,  delta: float, tau: float) -> float: ...

@dataclass
class HelmholtzParams:
    # Parâmetros do fluido
    M: float          # kg/mol
    Tc: float         # K
    rhoc_molar: float # mol/m^3 (atenção: base molar p/ δ=ρ/ρc)
    Rm: float = R_u   # J/(mol.K)

@dataclass
class Alpha0Ideal:
    """
    Parte ideal α0(δ,τ). Na forma short, α0 independe de δ (exceto ln δ em alguns modelos);
    aqui implementamos a forma típica: α0(τ) + ln(δ) (quando aplicável).
    """
    has_ln_delta: bool
    # componentes de α0(τ) dados por combinação de termos log/polinômios/exponenciais:
    # aqui usamos uma forma flexível via função de avaliação injetada
    def eval_all(self, delta: float, tau: float) -> Tuple[float,float,float]:
        """
        Retorna: (α0, α0_δ, α0_τ)
        Observação: no short-form polar usado por Lemmon & Span (2006), α0_δ = 1/δ (quando houver ln δ)
        """
        raise NotImplementedError("Fornecer implementação concreta para o fluido.")

class HelmholtzEOS:
    """
    EOS fundamental explícita em Helmholtz:
      φ(δ,τ) = α(δ,τ) = α0(δ,τ) + αr(δ,τ)
    com δ = ρ/ρc, τ = Tc/T ; ρ é densidade molar.
    """

    def __init__(self, params: HelmholtzParams, alpha0: Alpha0Ideal, residual_terms: list[ResidualTerm]):
        self.p = params
        self.alpha0 = alpha0
        self.terms = residual_terms

    # ------- blocos auxiliares (somatório residual) -------
    def _alphar_all(self, delta: float, tau: float):
        ar   = 0.0
        ar_d = 0.0
        ar_t = 0.0
        ar_dd= 0.0
        ar_dt= 0.0
        ar_tt= 0.0
        for term in self.terms:
            ar    += term.val(delta, tau)
            ar_d  += term.d_delta(delta, tau)
            ar_t  += term.d_tau(delta, tau)
            ar_dd += term.d_delta2(delta, tau)
            ar_dt += term.d_deltatau(delta, tau)
            ar_tt += term.d_tau2(delta, tau)
        return ar, ar_d, ar_t, ar_dd, ar_dt, ar_tt

    # =======================
    # Propriedades principais
    # =======================

    def pressure(self, rho_mass: float, T: float) -> float:
        """
        p [Pa] a partir de ρ_mass [kg/m3], T [K].
        Equações-padrão do formalismo Helmholtz (base molar → converter).
        """
        # converter ρ_mass -> ρ_molar
        rho_m = rho_mass / self.p.M
        delta = rho_m / self.p.rhoc_molar
        tau   = self.p.Tc / T

        # α0 e αr
        a0, a0_d, a0_t = self.alpha0.eval_all(delta, tau)
        ar, ar_d, ar_t, *_ = self._alphar_all(delta, tau)

        # p/(ρ R T) = 1 + δ*αr_δ  (forma padrão; α0 não contribui em δ)
        Z = 1.0 + delta * ar_d
        return rho_m * self.p.Rm * T * Z

    def enthalpy_mass(self, rho_mass: float, T: float) -> float:
        """
        h [J/kg] = (R_m T)*(1 + δ αr_δ + τ (α0_τ + αr_τ)) / M
        """
        rho_m = rho_mass / self.p.M
        delta = rho_m / self.p.rhoc_molar
        tau   = self.p.Tc / T

        a0, a0_d, a0_t = self.alpha0.eval_all(delta, tau)
        ar, ar_d, ar_t, *_ = self._alphar_all(delta, tau)

        h_molar = self.p.Rm * T * (1.0 + delta*ar_d + tau*(a0_t + ar_t))
        return h_molar / self.p.M