# params_N2O_lemmon_span_2006.py
from __future__ import annotations
from typing import List, Tuple
from dataclasses import dataclass
from eos_helmholtz_core import HelmholtzParams, HelmholtzEOS
from eos_helmholtz_short_forms import Alpha0_ShortPolar, PolyTerm, ExpoTerm
import math

# -----------------------------
#  Parâmetros críticos e M (N2O)
#  (valores de referência)
# -----------------------------
Tc   = 309.52067823146285        # K
pc   = 7.24481670160716e6        # Pa (não usado diretamente aqui)
rhoc_molar = 10290.903501226569  # mol/m^3
M    = 44.0128e-3                # kg/mol (0.0440128)
# Fonte (CoolProp doc / dados do artigo) [2](https://coolprop.org/fluid_properties/fluids/NitrousOxide.html)

params = HelmholtzParams(M=M, Tc=Tc, rhoc_molar=rhoc_molar)

# -----------------------------
#  α0(τ): parte ideal (short/polar)
#  Forma típica: α0(τ) = a + b*τ + ∑ c_i ln(1 - exp(-θ_i τ))
#  + ln(δ) (opcional)
# -----------------------------

# Estes coeficientes (a, b, {c_i, θ_i}) são da forma "short" polar (Span & Wagner; Lemmon & Span).
# A seguir, um exemplo de estrutura; os números exatos vêm do paper.
# >>> COLE AQUI (a, b, c_i, theta_i) conforme a TABELA do N2O (Lemmon & Span, 2006) <<<
a0_a = 0.0
a0_b = 0.0
a0_c = []       # lista de c_i
a0_th = []      # lista de θ_i

def a0_eval(tau: float) -> Tuple[float,float]:
    """
    Retorna α0(τ) e α0_τ(τ).
    α0(τ) = a + b τ + Σ c_i ln(1 - exp(-θ_i τ))
    α0_τ  = b + Σ c_i * [ ( -θ_i * exp(-θ_i τ) ) / (1 - exp(-θ_i τ)) ]
    """
    s = a0_a + a0_b * tau
    ds = a0_b
    for ci, thetai in zip(a0_c, a0_th):
        e = math.exp(-thetai*tau)
        s  += ci * math.log(1.0 - e)
        ds += ci * ( (-thetai*e) / (1.0 - e) )
    return s, ds

alpha0 = Alpha0_ShortPolar(has_ln_delta=True, a0_eval_func=a0_eval)

# -----------------------------
#  αr(δ,τ): parte residual "short" polar
#  Conjunto de termos (polinomiais + exponenciais)
# -----------------------------

# 1) EXPOENTES FIXOS da forma polar (definidos pela "classe" do modelo short)
# >>> COLE AQUI (listas d_i, t_i, c_i, alpha_i, beta_i, gamma_i, epsilon_i) conforme o paper <<<

# 2) COEFICIENTES n_i do N2O (12 termos no total)
# >>> COLE AQUI (lista n_i) conforme a TABELA do N2O (Lemmon & Span, 2006) <<<

# Abaixo, demonstro a montagem; substitua os arrays [] pelos valores do paper.
n_list      : List[float] = [0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, -0.0036985, -0.23263, -0.00042859, -0.042810]          # 12 coeficientes
d_list      : List[float] = []          # expoentes δ
t_list      : List[float] = []          # expoentes τ
c_list      : List[float] = []          # pode ser 0 na short polar
alpha_list  : List[float] = []          # parâmetros dos termos exponenciais
beta_list   : List[float] = []
gamma_list  : List[float] = []
eps_list    : List[float] = []

def build_terms():
    terms = []
    # Por convenção, nas formas short polares, parte dos 12 termos é POLY e parte EXPO
    # A divisão (quais índices são poly vs expo) é dada pelo artigo.
    # Exemplo ilustrativo (AJUSTE conforme a tabela do N2O):
    #   primeiros Np termos -> PolyTerm ; restantes -> ExpoTerm
    Np = 0  # <<< AJUSTE: número de termos polinomiais
    for i in range(len(n_list)):
        if i < Np:
            terms.append(PolyTerm(n=n_list[i], d=d_list[i], t=t_list[i], c=c_list[i] if i < len(c_list) else 0.0))
        else:
            terms.append(ExpoTerm(n=n_list[i], d=d_list[i], t=t_list[i],
                                  alpha=alpha_list[i - Np], beta=beta_list[i - Np],
                                  gamma=gamma_list[i - Np], epsilon=eps_list[i - Np]))
    return terms

def make_eos() -> HelmholtzEOS:
    terms = build_terms()
    return HelmholtzEOS(params=params, alpha0=alpha0, residual_terms=terms)